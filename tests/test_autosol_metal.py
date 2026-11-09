import io
import json
import stat
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.autosol import AutoSolPreparationError
from nasolve.autosol_metal import (
    FAILED,
    FOUND,
    NOT_FOUND,
    execute_metal_search,
    parse_metals,
)
from nasolve.cli import build_parser, main
from nasolve.config import AppConfig
from nasolve.coot_runtime import CootInstallation
from nasolve.coot_view import CootViewError, launch_coot_view, resolve_view_profile
from nasolve.model_assessment import file_sha256
from nasolve.phenix_runtime import PhenixInstallation
from nasolve.run_context import resolve_artifact_path

from .helpers import pdb_record
from .test_autosol import make_mtz_dump as _autosol_mtz_dump


def make_mtz_dump(root: Path) -> Path:
    # Same beamline wavelength as the fixture summary.html (0.9795 A).
    return _autosol_mtz_dump(root, wavelength=0.9795)


ENV = {"PATH": "/usr/bin:/bin"}
# Official AutoSol PHIL paths accepted by the fake program. Any other key
# (a typo, a short ambiguous name, or "sites") makes it fail like Phenix.
REQUIRED_KEYS = {
    "autosol.data", "autosol.seq_file", "autosol.labels", "autosol.atom_type",
    "autosol.lambda", "autosol.phasing.input_partpdb_file",
    "autosol.phasing.partpdb_rms", "autosol.model_building.build",
    "autosol.model_building.phase_improve_and_build", "autosol.general.nproc",
}
OPTIONAL_KEYS = {"autosol.phasing.mad_ha_add_list"}


def cryst1(space_group: str = "P 1") -> str:
    return f"CRYST1{50.0:9.3f}{50.0:9.3f}{60.0:9.3f}{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} {space_group:<11s}   1\n"


def make_autosol(root: Path, sites=(), return_code: int = 0, eff_rmsd: str | None = None) -> Path:
    """Fake phenix.autosol writing ``overall_best_ha_pdb.pdb`` with ``sites``.

    ``sites`` is a sequence of (element, (x, y, z)); an empty sequence writes
    no heavy-atom model at all.
    """
    records = "".join(
        pdb_record("HETATM", index, element, element, "A", index, occupancy=0.71,
                   element=element, x=xyz[0], y=xyz[1], z=xyz[2])
        for index, (element, xyz) in enumerate(sites, start=1)
    )
    executable = root / "phenix.autosol"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, sys\n"
        "from pathlib import Path\n"
        f"required = {sorted(REQUIRED_KEYS)!r}\n"
        f"optional = {sorted(OPTIONAL_KEYS)!r}\n"
        "args = dict(arg.partition('=')[::2] for arg in sys.argv[1:])\n"
        "keys = [arg.partition('=')[0] for arg in sys.argv[1:]]\n"
        "unknown = set(keys) - set(required) - set(optional)\n"
        "if unknown or not set(required) <= set(keys) or len(keys) != len(set(keys)):\n"
        "    print('Sorry: unknown or missing AutoSol parameters: ' + ', '.join(sorted(unknown)))\n"
        "    raise SystemExit(1)\n"
        "Path('received_args.json').write_text(json.dumps(sys.argv[1:]))\n"
        f"if {return_code}:\n"
        "    print('Sorry: no heavy-atom solution')\n"
        f"    raise SystemExit({return_code})\n"
        "out = Path('AutoSol_run_1_')\n"
        "out.mkdir()\n"
        f"rmsd = {eff_rmsd!r} or args['autosol.phasing.partpdb_rms']\n"
        "(out / 'autosol.eff').write_text(\n"
        "    'sites = None\\nbuild = False\\nphase_improve_and_build = False\\n'\n"
        "    f\"nproc = {args['autosol.general.nproc']}\\npartpdb_rms = {rmsd}\\n\")\n"
        f"records = {records!r}\n"
        "if records:\n"
        "    (out / 'overall_best_ha_pdb.pdb').write_text(records + 'END\\n')\n"
        "print('AutoSol finished')\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_run(root: Path, *, tfz: float = 12.3, space_group: str = "P 1", stage: str = "phaser") -> Path:
    dataset = root / "dataset"
    frame = root / "MR_frames" / "5W6W"
    run = dataset / "AutoMR" / "run_001"
    phaser = run / "Phaser"
    frame.mkdir(parents=True)
    phaser.mkdir(parents=True)
    reflections = dataset / "staraniso_alldata-unique.mtz"
    reflections.write_bytes(b"mtz")
    summary = dataset / "summary.html"
    summary.write_text("wavelength [A] = 0.979500\n")
    frame_model = frame / "C_G.pdb"
    frame_model.write_text("END\n")
    (frame / "seq_base.txt").write_text("GAGC\n\nCTGC\n")
    model = phaser / "mr_solution.pdb"
    model.write_text(
        cryst1(space_group)
        + pdb_record("ATOM", 1, "N3", "DC", "A", 12, element="N", x=-2.1, y=0.0, z=0.0)
        + pdb_record("ATOM", 2, "C1'", "DC", "A", 12, element="C", x=-6.0, y=0.0, z=0.0)
        + pdb_record("ATOM", 3, "N3", "DC", "B", 4, element="N", x=2.2, y=0.0, z=0.0)
        + pdb_record("ATOM", 4, "N7", "DG", "B", 5, element="N", x=9.0, y=0.0, z=0.0)
        + "END\n"
    )
    phaser_map = phaser / "mr_solution.mtz"
    phaser_map.write_bytes(b"phaser map coefficients")
    status = "MR_SUCCESS" if tfz >= 8 else ("MR_REVIEW" if tfz >= 7 else "MR_FAILED")
    report = {
        "workflow": "automr",
        "stage": stage,
        "status": status if stage == "phaser" else "POSTMR_READY",
        "message": "fixture",
        "inputs": {
            "reflections": str(reflections.resolve()),
            "reflections_sha256": file_sha256(reflections),
            "summary": str(summary.resolve()),
            "model": str(frame_model.resolve()),
        },
        "execution": {
            "phaser": {
                "solution_pdb": str(model.resolve()),
                "solution_mtz": str(phaser_map.resolve()),
                "best_tfz": tfz,
            },
        },
    }
    (run / "report.json").write_text(json.dumps(report))
    return run


def search(run: Path, autosol: Path, root: Path, metals=("AG",), **kwargs):
    return execute_metal_search(
        run, autosol, make_mtz_dump(root), metals,
        environment=ENV, processor_count=8, **kwargs,
    )


class ParseMetalTests(unittest.TestCase):
    def test_metals_are_normalized_ordered_and_deduplicated(self):
        self.assertEqual(parse_metals(["ag", "Hg,ZN", "AG"]), ("Ag", "Hg", "Zn"))

    def test_non_metals_and_empty_input_are_rejected(self):
        for bad in (["I"], ["Se"], ["XX"], [" , "]):
            with self.subTest(bad=bad), self.assertRaises(AutoSolPreparationError):
                parse_metals(bad)


class MetalSearchTests(unittest.TestCase):
    def test_found_site_runs_mr_sad_with_phaser_model_and_rmsd_one(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)

            self.assertEqual(result.status, FOUND)
            self.assertEqual(result.attempt_directory, run / "AutoSolMetal" / "attempt_001")
            args = json.loads((result.attempt_directory / "received_args.json").read_text())
            self.assertIn("autosol.atom_type=Ag", args)
            self.assertIn("autosol.phasing.partpdb_rms=1", args)
            self.assertIn("autosol.model_building.build=False", args)
            self.assertIn("autosol.model_building.phase_improve_and_build=False", args)
            self.assertFalse(any(a.partition("=")[0].endswith(".sites") for a in args))
            model_arg = next(a for a in args if a.startswith("autosol.phasing.input_partpdb_file="))
            self.assertTrue(model_arg.endswith("/Phaser/mr_solution.pdb"))
            self.assertIn("autosol.lambda=0.9795", args)

            site = result.sites[0]
            self.assertEqual(site["element"], "AG")
            nearest = site["nearest_residues"]
            self.assertEqual((nearest[0]["chain"], nearest[0]["residue_number"]), ("A", "12"))
            self.assertAlmostEqual(nearest[0]["distance"], 2.1)
            self.assertEqual((nearest[1]["chain"], nearest[1]["residue_number"]), ("B", "4"))

    def test_report_is_additive_and_portable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            report = json.loads((run / "report.json").read_text())

            # The run's own stage and status are untouched.
            self.assertEqual((report["stage"], report["status"]), ("phaser", "MR_SUCCESS"))
            block = report["metal_autosol"]
            self.assertEqual(block["schema_version"], 1)
            self.assertEqual(block["latest"], "attempt_001")
            attempt = block["attempts"][0]
            self.assertEqual(attempt["status"], FOUND)
            self.assertFalse(attempt["model_changed"])
            self.assertEqual(attempt["effective_parameters"]["partpdb_rms"], 1.0)
            ha = attempt["outputs"]["heavy_atom_model"]
            self.assertEqual(ha["anchor"], "run")
            self.assertTrue(ha["relative_path"].endswith("overall_best_ha_pdb.pdb"))
            self.assertRegex(ha["sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(attempt["inputs"]["phaser_map"]["relative_path"], "Phaser/mr_solution.mtz")
            self.assertEqual(resolve_artifact_path(ha, run), (run / ha["relative_path"]).resolve())
            self.assertTrue((run / "AutoSolMetal" / "attempt_001" / "report.json").is_file())

    def test_several_metals_and_custom_rmsd(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(
                run, make_autosol(root, [("HG", (0.0, 0.0, 0.0))]), root,
                metals=("AG", "hg"), partial_model_rmsd=0.8,
            )
            args = json.loads((result.attempt_directory / "received_args.json").read_text())
            self.assertIn("autosol.atom_type=Ag", args)
            self.assertIn("autosol.phasing.mad_ha_add_list=Hg", args)
            self.assertIn("autosol.phasing.partpdb_rms=0.8", args)
            self.assertEqual(result.status, FOUND)
            self.assertEqual(result.metals, ("Ag", "Hg"))

    def test_no_site_is_a_normal_not_found_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(run, make_autosol(root, []), root)
            self.assertEqual(result.status, NOT_FOUND)
            self.assertIn("normal", result.message)
            self.assertIsNone(result.heavy_atom_model)

    def test_only_other_elements_is_not_found(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(run, make_autosol(root, [("ZN", (0.0, 0.0, 0.0))]), root)
            self.assertEqual(result.status, NOT_FOUND)
            attempt = json.loads((run / "report.json").read_text())["metal_autosol"]["attempts"][0]
            self.assertEqual(attempt["other_sites"], 1)

    def test_autosol_error_is_recorded_not_raised(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(run, make_autosol(root, return_code=1), root)
            self.assertEqual(result.status, FAILED)
            self.assertIn("exit status 1", result.message)
            self.assertIn("no heavy-atom solution", result.log_path.read_text())

    def test_rmsd_not_applied_by_autosol_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(
                run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))], eff_rmsd="1.5"), root
            )
            self.assertEqual(result.status, FAILED)
            self.assertIn("partpdb_rms is 1.5", result.message)

    def test_mr_review_needs_permission_and_failed_mr_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, tfz=7.4)
            autosol = make_autosol(root, [("AG", (0.0, 0.0, 0.0))])
            with self.assertRaisesRegex(AutoSolPreparationError, "allow-mr-review"):
                search(run, autosol, root)
            self.assertFalse((run / "AutoSolMetal").exists())  # refused: no trace
            self.assertEqual(search(run, autosol, root, allow_mr_review=True).status, FOUND)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, tfz=5.0)
            with self.assertRaisesRegex(AutoSolPreparationError, "below 7.0"):
                search(run, make_autosol(root), root)

    def test_changed_reflections_stop_before_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            (root / "dataset" / "staraniso_alldata-unique.mtz").write_bytes(b"changed")
            with self.assertRaises(AutoSolPreparationError):
                search(run, make_autosol(root), root)
            self.assertFalse((run / "AutoSolMetal").exists())

    def test_repeated_searches_never_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            first = search(run, make_autosol(root, []), root)
            second = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            self.assertEqual(first.attempt_directory.name, "attempt_001")
            self.assertEqual(second.attempt_directory.name, "attempt_002")
            self.assertTrue((first.attempt_directory / "autosol.console.log").is_file())
            block = json.loads((run / "report.json").read_text())["metal_autosol"]
            self.assertEqual([a["status"] for a in block["attempts"]], [NOT_FOUND, FOUND])
            self.assertEqual(block["latest"], "attempt_002")

    def test_works_after_postmr_and_with_unconfigured_symmetry(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, stage="postmr", space_group="P 63")
            result = search(run, make_autosol(root, [("HG", (0.0, 0.0, 0.0))]), root, metals=("HG",))
            self.assertEqual(result.status, FOUND)
            attempt = json.loads((run / "report.json").read_text())["metal_autosol"]["attempts"][0]
            self.assertIn("same asymmetric unit only", attempt["distance_method"])
            self.assertEqual(json.loads((run / "report.json").read_text())["status"], "POSTMR_READY")


class MetalViewTests(unittest.TestCase):
    def test_metal_view_opens_phaser_model_sites_and_phaser_map(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            with self.assertRaisesRegex(CootViewError, "No metal search"):
                resolve_view_profile(run, stage="metal")
            result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)

            profile = resolve_view_profile(run, stage="metal")
            self.assertEqual(profile.model_path, (run / "Phaser" / "mr_solution.pdb").resolve())
            self.assertEqual(profile.map_path, (run / "Phaser" / "mr_solution.mtz").resolve())
            self.assertEqual(profile.extra_model_paths, (result.heavy_atom_model.resolve(),))

            coot = root / "coot"
            coot.write_text("#!/bin/sh\n")
            view = launch_coot_view(
                run, coot, stage="metal", environment=ENV,
                launcher=lambda *args, **kwargs: SimpleNamespace(pid=4321),
            )
            command = list(view.command)
            self.assertEqual(command[command.index("--auto") + 1], str(profile.map_path))
            pdbs = [command[i + 1] for i, value in enumerate(command) if value == "--pdb"]
            self.assertEqual(pdbs, [str(profile.model_path), str(result.heavy_atom_model.resolve())])
            self.assertEqual(view.working_directory, run / "CootGUI" / "metal")

    def test_view_refuses_a_changed_heavy_atom_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            result.heavy_atom_model.write_text("tampered\n")
            with self.assertRaises(CootViewError):
                resolve_view_profile(run, stage="metal")


class MetalCommandTests(unittest.TestCase):
    def _phenix(self, root: Path, autosol: Path) -> PhenixInstallation:
        return PhenixInstallation(
            root=None, setup_script=None, version="2.2.1-6174",
            executables={"phenix.autosol": autosol, "phenix.mtz.dump": make_mtz_dump(root)},
            environment=ENV, source="test",
        )

    def _main(self, argv, phenix, **patches):
        output, errors = io.StringIO(), io.StringIO()
        with (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli.discover_phenix", return_value=phenix),
            patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
            patch("nasolve.cli.remember_coot"),
            patch("nasolve.cli.discover_coot", return_value=patches.get("coot")),
            patch("nasolve.cli.launch_coot_view", side_effect=patches.get("launch")) as launch,
            redirect_stdout(output), redirect_stderr(errors),
        ):
            code = main(argv)
        return code, output.getvalue(), errors.getvalue(), launch

    def test_parser_accepts_repeated_metals(self):
        args = build_parser().parse_args(["autosol", "run_001", "--metal", "AG", "--metal", "hg"])
        self.assertEqual(args.metal, ["AG", "hg"])
        self.assertIsNone(args.rmsd)
        self.assertEqual(build_parser().parse_args(["show", "run", "--stage", "metal"]).stage, "metal")

    def test_found_sites_open_coot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            coot = root / "coot"
            coot.write_text("#!/bin/sh\n")
            launch = partial(launch_coot_view, launcher=lambda *a, **k: SimpleNamespace(pid=4321))
            code, out, _, launched = self._main(
                ["autosol", str(run), "--metal", "AG"],
                self._phenix(root, make_autosol(root, [("AG", (0.0, 0.0, 0.0))])),
                coot=CootInstallation(coot, "1.3.3", "test"), launch=launch,
            )
            self.assertEqual(code, 0, out)
            self.assertIn(FOUND, out)
            self.assertIn("A:12 DC N3 2.10 A", out)
            self.assertIn("Opened Coot (PID 4321)", out)
            self.assertEqual(launched.call_args.kwargs["stage"], "metal")
            self.assertTrue((run / "CootGUI" / "metal" / "launch.json").is_file())

    def test_not_found_and_no_coot_do_not_open_coot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            code, out, _, launched = self._main(
                ["autosol", str(run), "--metal", "AG"], self._phenix(root, make_autosol(root, []))
            )
            self.assertEqual(code, 0)
            self.assertIn(NOT_FOUND, out)
            self.assertIn("nothing was added to the model", out)
            launched.assert_not_called()

            code, out, _, launched = self._main(
                ["autosol", str(run), "--metal", "AG", "--no-coot"],
                self._phenix(root, make_autosol(root, [("AG", (0.0, 0.0, 0.0))])),
            )
            self.assertEqual(code, 0)
            self.assertIn("--stage metal", out)
            launched.assert_not_called()

    def test_metal_options_without_metal_and_bad_metal_are_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            phenix = self._phenix(root, make_autosol(root))
            code, _, errors, _ = self._main(["autosol", str(run), "--rmsd", "0.8"], phenix)
            self.assertEqual(code, 2)
            self.assertIn("add --metal", errors)
            code, _, errors, _ = self._main(["autosol", str(run), "--metal", "I"], phenix)
            self.assertEqual(code, 2)
            self.assertIn("not a metal", errors)
            self.assertFalse((run / "AutoSolMetal").exists())


if __name__ == "__main__":
    unittest.main()
