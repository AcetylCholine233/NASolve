import io
import json
import os
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
    model_sequences,
    FAILED,
    FOUND,
    NOT_FOUND,
    execute_metal_search,
    parse_metals,
)
from nasolve.cli import build_parser, main
from nasolve.config import AppConfig
from nasolve.coot_runtime import CootDiscoveryError, CootInstallation
from nasolve.coot_view import CootViewError, launch_coot_view, resolve_run, resolve_view_profile
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
            # Regression: searches were written inside the AutoMR run.
            dataset = root / "dataset"
            self.assertEqual(result.output_directory, dataset / "AutoSolMetal" / "run_001")
            self.assertFalse((run / "AutoSolMetal").exists())
            self.assertEqual(result.run_directory, run.resolve())
            args = json.loads((result.output_directory / "received_args.json").read_text())
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
            result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            report = json.loads((run / "report.json").read_text())
            metal = result.output_directory

            # The run's own stage and status are untouched.
            self.assertEqual((report["stage"], report["status"]), ("phaser", "MR_SUCCESS"))
            block = report["metal_autosol"]
            self.assertEqual(block["schema_version"], 2)
            self.assertEqual(block["latest"], "run_001")
            listed = block["searches"][0]
            self.assertEqual(listed["status"], FOUND)
            self.assertEqual(listed["directory"],
                             {"anchor": "dataset", "relative_path": "AutoSolMetal/run_001"})
            self.assertEqual(resolve_artifact_path(listed["report"], run), metal / "report.json")

            # The full record lives with the search, anchored to it or the dataset.
            record = json.loads((metal / "report.json").read_text())
            self.assertEqual((record["workflow"], record["metal_run"]), ("metal_autosol", "run_001"))
            self.assertEqual(record["source_run"],
                             {"anchor": "dataset", "relative_path": "AutoMR/run_001"})
            self.assertFalse(record["model_changed"])
            self.assertEqual(record["effective_parameters"]["partpdb_rms"], 1.0)
            ha = record["outputs"]["heavy_atom_model"]
            self.assertEqual(ha["anchor"], "run")
            self.assertTrue(ha["relative_path"].endswith("overall_best_ha_pdb.pdb"))
            self.assertRegex(ha["sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(record["inputs"]["phaser_map"],
                             {**record["inputs"]["phaser_map"], "anchor": "dataset",
                              "relative_path": "AutoMR/run_001/Phaser/mr_solution.mtz"})
            self.assertEqual(resolve_artifact_path(ha, metal), (metal / ha["relative_path"]).resolve())

    def test_metal_search_record_survives_moving_the_dataset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            moved = root / "moved" / "dataset"
            moved.parent.mkdir()
            (root / "dataset").rename(moved)
            profile = resolve_view_profile(moved / "AutoSolMetal" / "run_001")
            self.assertEqual(profile.model_path,
                             (moved / "AutoMR" / "run_001" / "Phaser" / "mr_solution.pdb").resolve())
            self.assertEqual(
                resolve_view_profile(moved / "AutoMR" / "run_001", stage="metal").extra_model_paths,
                profile.extra_model_paths,
            )

    def test_several_metals_and_custom_rmsd(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(
                run, make_autosol(root, [("HG", (0.0, 0.0, 0.0))]), root,
                metals=("AG", "hg"), partial_model_rmsd=0.8,
            )
            args = json.loads((result.output_directory / "received_args.json").read_text())
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
            record = json.loads(result.report_path.read_text())
            self.assertEqual(record["other_sites"], 1)

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
            self.assertFalse((root / "dataset" / "AutoSolMetal").exists())  # refused: no trace
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
            self.assertFalse((root / "dataset" / "AutoSolMetal").exists())

    def test_repeated_searches_never_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            first = search(run, make_autosol(root, []), root)
            second = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            self.assertEqual(first.output_directory.name, "run_001")
            self.assertEqual(second.output_directory.name, "run_002")
            self.assertTrue((first.output_directory / "autosol.console.log").is_file())
            block = json.loads((run / "report.json").read_text())["metal_autosol"]
            self.assertEqual([a["status"] for a in block["searches"]], [NOT_FOUND, FOUND])
            self.assertEqual(block["latest"], "run_002")

    def test_works_after_postmr_and_with_unconfigured_symmetry(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, stage="postmr", space_group="P 63")
            result = search(run, make_autosol(root, [("HG", (0.0, 0.0, 0.0))]), root, metals=("HG",))
            self.assertEqual(result.status, FOUND)
            record = json.loads(result.report_path.read_text())
            self.assertIn("same asymmetric unit only", record["distance_method"])
            self.assertEqual(json.loads((run / "report.json").read_text())["status"], "POSTMR_READY")


class MetalViewTests(unittest.TestCase):
    def test_metal_view_opens_phaser_model_sites_and_phaser_map(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            with self.assertRaisesRegex(CootViewError, "No metal search"):
                resolve_view_profile(run, stage="metal")
            not_found = search(run, make_autosol(root, []), root)
            with self.assertRaisesRegex(CootViewError, "did not find sites"):
                resolve_view_profile(not_found.output_directory)
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
            # Opening the metal-search folder itself gives the same view.
            direct = resolve_view_profile(result.output_directory)
            self.assertEqual(
                (direct.model_path, direct.map_path, direct.extra_model_paths),
                (profile.model_path, profile.map_path, profile.extra_model_paths),
            )
            with self.assertRaisesRegex(CootViewError, "only be shown as the metal view"):
                resolve_view_profile(result.output_directory, stage="automr")

    def test_schema_1_attempt_inside_the_run_is_still_viewable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            attempt = run / "AutoSolMetal" / "attempt_001"
            attempt.mkdir(parents=True)
            sites = attempt / "overall_best_ha_pdb.pdb"
            sites.write_text("HETATM    1 AG    AG A   1       0.000   0.000   0.000  1.00 20.00          AG\n")
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["metal_autosol"] = {"schema_version": 1, "latest": "attempt_001", "attempts": [{
                "attempt": "attempt_001", "status": FOUND,
                "inputs": {
                    "phaser_model": {"anchor": "run", "relative_path": "Phaser/mr_solution.pdb"},
                    "phaser_map": {"anchor": "run", "relative_path": "Phaser/mr_solution.mtz"},
                },
                "outputs": {"heavy_atom_model": {
                    "anchor": "run", "relative_path": "AutoSolMetal/attempt_001/overall_best_ha_pdb.pdb"}},
            }]}
            report_path.write_text(json.dumps(report))
            profile = resolve_view_profile(run, stage="metal")
            self.assertEqual(profile.extra_model_paths, (sites.resolve(),))
            self.assertIn("attempt_001", profile.source)

    def test_show_finds_a_metal_search_from_the_dataset_or_the_active_dataset(self):
        # Regression: "show AutoSolMetal/run_001" from outside the dataset said
        # "AutoMR report is missing" without explaining what was expected.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            run = make_run(root)
            result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            dataset = root / "dataset"
            elsewhere = root / "elsewhere"
            elsewhere.mkdir()
            previous = os.getcwd()
            try:
                os.chdir(dataset)
                self.assertEqual(resolve_run("AutoSolMetal/run_001"), result.output_directory)
                os.chdir(elsewhere)
                with self.assertRaisesRegex(
                    CootViewError, r"No NASolve report\.json in .*elsewhere.*AutoSolMetal/run_NNN"
                ):
                    resolve_run("AutoSolMetal/run_001")
                self.assertEqual(
                    resolve_run("AutoSolMetal/run_001", active_dataset=dataset),
                    result.output_directory,
                )

                config = AppConfig()
                config.workspace.dataset = str(dataset)
                output = io.StringIO()
                with (
                    patch("nasolve.cli.load_config", return_value=config),
                    patch("nasolve.cli.save_config"), patch("nasolve.cli.remember_coot"),
                    patch("nasolve.cli.discover_coot",
                          return_value=CootInstallation(root / "coot", "1.3.3", "test")),
                    patch("nasolve.cli.launch_coot_view", side_effect=partial(
                        launch_coot_view, launcher=lambda *a, **k: SimpleNamespace(pid=7))),
                    redirect_stdout(output), redirect_stderr(io.StringIO()),
                ):
                    (root / "coot").write_text("#!/bin/sh\n")
                    code = main(["show", "AutoSolMetal/run_001"])
                self.assertEqual(code, 0, output.getvalue())
                self.assertIn("Opened metal in run: " + str(result.output_directory), output.getvalue())
                self.assertIn("Merge metal sites and save", output.getvalue())
            finally:
                os.chdir(previous)

    def test_view_refuses_a_changed_heavy_atom_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root)
            result.heavy_atom_model.write_text("tampered\n")
            with self.assertRaises(CootViewError):
                resolve_view_profile(run, stage="metal")


class SequenceSourceTests(unittest.TestCase):
    """AutoSol's sequence for the metal search: --sequence, the run, or the
    Phaser model. Regression: a run without MR_frames/FRAME/seq_base.txt (the
    3GBI frame has none) stopped with "AutoSol needs ... seq_base.txt"."""

    def test_sequence_from_the_phaser_model_when_the_run_has_none(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            (root / "MR_frames" / "5W6W" / "seq_base.txt").unlink()
            result = search(run, make_autosol(root, []), root)
            record = json.loads(result.report_path.read_text())
            # Only A:12 has a C1' atom in this fixture; B:4 and B:5 are base
            # fragments, so they are not counted as nucleotides.
            self.assertEqual(record["sequence_source"], "derived from the Phaser model (chain A: 1)")
            self.assertEqual(record["sequence_warnings"], [])
            sequence = result.output_directory / "sequence_input.txt"
            self.assertEqual(sequence.read_text(), "C\n")
            self.assertIn(f"autosol.seq_file={sequence}", record["command"])

    def test_run_sequence_is_still_preferred_and_sequence_file_overrides(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            result = search(run, make_autosol(root, []), root)
            record = json.loads(result.report_path.read_text())
            self.assertEqual(record["sequence_source"], "run-frozen sequence file")
            self.assertEqual((result.output_directory / "sequence_input.txt").read_text(),
                             "GAGC\n\nCTGC\n")

            given = root / "my_sequence.txt"
            given.write_text(">strand 1\nGAGCAG\n")
            result = search(run, make_autosol(root, []), root, sequence_file=given)
            record = json.loads(result.report_path.read_text())
            self.assertEqual(record["sequence_source"], "--sequence my_sequence.txt")
            self.assertEqual((result.output_directory / "sequence_input.txt").read_text(),
                             ">strand 1\nGAGCAG\n")

            (root / "empty.txt").write_text(">only a header\n\n")
            for path, message in ((root / "missing.txt", "No sequence file"),
                                  (root / "empty.txt", "has no sequence")):
                with self.subTest(path=path.name), self.assertRaisesRegex(
                    AutoSolPreparationError, message
                ):
                    search(run, make_autosol(root, []), root, sequence_file=path)

    def test_no_sequence_anywhere_names_the_option(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            (root / "MR_frames" / "5W6W" / "seq_base.txt").unlink()
            (run / "Phaser" / "mr_solution.pdb").write_text(
                pdb_record("HETATM", 1, "AG", "AG", "B", 101, element="AG") + "END\n"
            )
            with self.assertRaisesRegex(AutoSolPreparationError, "give one with --sequence FILE"):
                search(run, make_autosol(root, []), root)
            self.assertFalse((root / "dataset" / "AutoSolMetal").exists())

    def test_model_sequences(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.pdb"
            records = []
            for serial, (name, chain, number) in enumerate((
                ("DG", "A", 1), ("DU", "A", 2), ("5IU", "A", 3), ("XYZ", "A", 4),
                ("DC", "B", 1), ("G", "R", 1), ("U", "R", 2),
            ), start=1):
                records.append(pdb_record("ATOM", serial, "C1'", name, chain, number, element="C"))
            records.append(pdb_record("HETATM", 20, "AG", "AG", "B", 101, element="AG"))
            records.append(pdb_record("HETATM", 21, "O", "HOH", "W", 1, element="O"))
            path.write_text("".join(records) + "END\n")
            chains, warnings = model_sequences(path)
            self.assertEqual(chains, {"A": "GTT", "B": "C", "R": "GU"})
            self.assertEqual(len(warnings), 1)
            self.assertIn("A:4 XYZ left out", warnings[0])


def with_pair(run: Path, codes=("C38", "5IU"), requested="iC:iU") -> None:
    """Give the AutoMR fixture a W-frame standard pair, as AutoMR records it."""
    report_path = run / "report.json"
    report = json.loads(report_path.read_text())
    report["frame"] = {"name": "W", "catalogue_directory": "5W6W"}
    report["post_mr_plan"] = {
        "standard_pair": {"requested": requested, "ligand_codes": list(codes)},
        "mutations": {}, "sequences": {},
    }
    report_path.write_text(json.dumps(report))


def fake_build(original, plan, model_dir, coot_dir, coot_executable, environment, *,
               output_name="after_coot.pdb"):
    """Stand-in for postmr.apply_model_mutations: writes the planned residues,
    with the iodine of each 5-iodo base where AutoSol found the side sites."""
    calls.append({"plan": plan, "coot": coot_executable, "output_name": output_name,
                  "original": original})
    lines = []
    for serial, action in enumerate(plan.coot_actions, start=1):
        chain, number = action.site.split(":")
        x = -2.1 if chain == "A" else 2.2
        lines.append(pdb_record("HETATM", serial, "C1'", action.after, chain, int(number),
                                element="C", x=x * 3, y=0.0, z=0.0))
        lines.append(pdb_record("HETATM", serial + 10, "I", action.after, chain, int(number),
                                element="I", x=x, y=4.0, z=0.0))
    built = model_dir / output_name
    built.write_text(cryst1("P 1") + "".join(lines) + "END\n")
    log = coot_dir / "coot.log"
    log.write_text("NASOLVE_RESULT\n")
    return built, log, {"restored": [action.site for action in plan.coot_actions]}


calls: list[dict] = []


class BuildRequestedResiduesTests(unittest.TestCase):
    """Regression (iCiU_Ag): the metal search used the plain Phaser model, so
    the iodines of 5-iodo-dC and 5-iodo-dU were found as two extra "Ag" sites.
    The AutoMR-requested residues are now built in first, as in PostMR."""

    def setUp(self):
        calls.clear()

    def test_requested_pair_is_built_into_the_partial_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            with_pair(run)
            coot = root / "coot"
            coot.write_text("#!/bin/sh\n")
            with patch("nasolve.postmr.apply_model_mutations", side_effect=fake_build):
                result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root,
                                coot_executable=coot)
            self.assertEqual(result.status, FOUND)
            (call,) = calls
            self.assertEqual(call["coot"], coot)
            self.assertEqual(call["output_name"], "partial_model.pdb")
            self.assertEqual(
                [(a.site, a.before, a.after, a.method) for a in call["plan"].coot_actions],
                [("A:12", "DC", "C38", "coot-parent-overlap"),
                 ("B:4", "DC", "5IU", "coot-parent-overlap")],
            )
            partial = result.output_directory / "PartialModel" / "partial_model.pdb"
            self.assertEqual(result.partial_model, partial)
            self.assertEqual(call["original"].read_bytes(), result.phaser_model.read_bytes())
            record = json.loads(result.report_path.read_text())
            self.assertIn(f"autosol.phasing.input_partpdb_file={partial}", record["command"])
            self.assertEqual(record["partial_model"]["mutations"], [
                {"site": "A:12", "before": "DC", "after": "C38", "method": "coot-parent-overlap"},
                {"site": "B:4", "before": "DC", "after": "5IU", "method": "coot-parent-overlap"},
            ])
            self.assertIn("AutoMR residues built in", record["partial_model"]["source"])
            self.assertEqual(record["partial_model"]["coordinate_restoration"],
                             {"restored": ["A:12", "B:4"]})
            self.assertEqual(record["inputs"]["partial_model"]["relative_path"],
                             "PartialModel/partial_model.pdb")
            self.assertFalse(record["model_changed"])
            # Distances are measured against the built model.
            nearest = record["sites"][0]["nearest_residues"][0]
            self.assertIn(nearest["residue_name"], {"C38", "5IU"})
            # The Coot view shows the built model; the merge keeps the Phaser name.
            profile = resolve_view_profile(result.output_directory)
            self.assertEqual(profile.model_path, partial.resolve())
            self.assertEqual(profile.metal_merge.merged_stem, "mr_solution-merged")
            self.assertIn("AutoMR residues built in", profile.source)

    def test_coot_is_needed_only_when_there_is_something_to_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            with_pair(run)
            with self.assertRaisesRegex(AutoSolPreparationError,
                                        r"asks for A:12 DC->C38, B:4 DC->5IU.*no Coot"):
                search(run, make_autosol(root, []), root)
            self.assertFalse((root / "dataset" / "AutoSolMetal").exists())
            # --original-model: the Phaser model as it is, no Coot.
            result = search(run, make_autosol(root, []), root, build_requested_residues=False)
            self.assertEqual(result.partial_model, result.phaser_model)
            record = json.loads(result.report_path.read_text())
            self.assertEqual(record["partial_model"]["source"],
                             "Phaser model as it is (--original-model)")
            self.assertEqual(calls, [])

    def test_a_pair_that_cannot_be_built_is_named(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            # "IC:IU" resolves to isocytidine and 5-iodouridine (RNA), not the
            # curated 5-iodo-dC/dU (written iC:iU).
            with_pair(run, codes=("IC", "IU"), requested="IC:IU")
            with self.assertRaisesRegex(AutoSolPreparationError,
                                        r"Cannot build the AutoMR residues.*IC.*--original-model"):
                search(run, make_autosol(root, []), root, coot_executable=root / "coot")
            self.assertFalse((root / "dataset" / "AutoSolMetal").exists())

    def test_coot_failure_stops_before_autosol(self):
        from nasolve.postmr import PostMRPreparationError

        def broken(*args, **kwargs):
            raise PostMRPreparationError("Coot mutation failed with status 1; inspect coot.log")

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            with_pair(run)
            with patch("nasolve.postmr.apply_model_mutations", side_effect=broken):
                result = search(run, make_autosol(root, [("AG", (0.0, 0.0, 0.0))]), root,
                                coot_executable=root / "coot")
            self.assertEqual(result.status, FAILED)
            self.assertIn("Coot could not build the AutoMR residues", result.message)
            self.assertFalse(result.log_path.exists())   # AutoSol never ran


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
            patch("nasolve.cli.discover_coot",
                  return_value=patches.get("coot"),
                  side_effect=None if patches.get("coot") else CootDiscoveryError("no Coot here")),
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
            metal = root / "dataset" / "AutoSolMetal" / "run_001"
            self.assertIn(f"Metal-search folder: {metal}", out)
            self.assertTrue((metal / "CootGUI" / "metal" / "launch.json").is_file())
            self.assertFalse((run / "CootGUI").exists())

    def test_sequence_option(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            (root / "MR_frames" / "5W6W" / "seq_base.txt").unlink()
            code, out, err, _ = self._main(
                ["autosol", str(run), "--metal", "AG", "--no-coot"],
                self._phenix(root, make_autosol(root, [])),
            )
            self.assertEqual(code, 0, err)
            self.assertIn("Sequence: derived from the Phaser model (chain A: 1)", out)
            given = root / "seq.txt"
            given.write_text("GAGC\n")
            code, out, err, _ = self._main(
                ["autosol", str(run), "--metal", "AG", "--no-coot", "--sequence", str(given)],
                self._phenix(root, make_autosol(root, [])),
            )
            self.assertEqual(code, 0, err)
            self.assertIn("Sequence: --sequence seq.txt", out)
            code, _, err, _ = self._main(
                ["autosol", str(run), "--sequence", str(given)],
                self._phenix(root, make_autosol(root, [])),
            )
            self.assertEqual(code, 2)
            self.assertIn("--sequence, --original-model and --no-coot apply only to the "
                          "metal search", err)

    def test_original_model_option(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            with_pair(run)
            code, out, err, _ = self._main(
                ["autosol", str(run), "--metal", "AG", "--no-coot"],
                self._phenix(root, make_autosol(root, [])),
            )
            self.assertEqual(code, 2)
            self.assertIn("no Coot is configured", err)
            code, out, err, _ = self._main(
                ["autosol", str(run), "--metal", "AG", "--no-coot", "--original-model"],
                self._phenix(root, make_autosol(root, [])),
            )
            self.assertEqual(code, 0, err)
            self.assertIn("Partial model: Phaser model as it is (--original-model)", out)
            coot = root / "coot"
            coot.write_text("#!/bin/sh\n")
            with patch("nasolve.postmr.apply_model_mutations", side_effect=fake_build):
                code, out, err, _ = self._main(
                    ["autosol", str(run), "--metal", "AG", "--no-coot"],
                    self._phenix(root, make_autosol(root, [])),
                    coot=CootInstallation(coot, "1.3.3", "test"),
                )
            self.assertEqual(code, 0, err)
            self.assertIn("Partial model: Phaser model with the AutoMR residues built in", out)
            self.assertIn("  A:12 DC -> C38", out)
            self.assertIn("  B:4 DC -> 5IU", out)

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
            self.assertIn("Open later with: nasolve show", out)
            self.assertIn("AutoSolMetal", out)
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
            self.assertFalse((root / "dataset" / "AutoSolMetal").exists())


if __name__ == "__main__":
    unittest.main()
