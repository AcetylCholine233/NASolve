import json
import math
import shutil
import stat
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from nasolve.autosol import (
    AutoSolPreparationError,
    _read_run_report,
    _sequence_source,
    execute_autosol,
    nearest_symmetry_distance,
    read_wavelength,
)
from nasolve.model_assessment import file_sha256
from nasolve.run_context import artifact_reference, resolve_artifact_path

from .helpers import pdb_record, symlinked_temporary_directory


class AutoSolReportTests(unittest.TestCase):
    def test_non_object_report_is_rejected_cleanly(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "run_001"
            run.mkdir()
            (run / "report.json").write_text(json.dumps([]))

            with self.assertRaisesRegex(AutoSolPreparationError, "not a JSON object"):
                _read_run_report(run)


def cryst1(
    a: float = 50.0,
    b: float = 50.0,
    c: float = 60.0,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 120.0,
    space_group: str = "R 3",
) -> str:
    return (
        f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}"
        f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} {space_group:<11s}   9\n"
    )


def hex_cart(x: float, y: float, z: float) -> tuple[float, float, float]:
    return (50.0 * x - 25.0 * y, 50.0 * math.sin(math.radians(120.0)) * y, 60.0 * z)


def make_mtz_dump(root: Path, wavelength: float = 1.377618) -> Path:
    executable = root / "phenix.mtz.dump"
    # Header layout of phenix.mtz.dump: HKL_base (placeholder wavelength)
    # then the data dataset that carries the beamline wavelength.
    header = (
        "Crystal 1:\n  Dataset 1:\n    Name: HKL_base\n    Wavelength: 9.999\n"
        "    label #valid %valid min max type\n"
        "    H 10 100.00% 0 9 H: index h,k,l\n"
        "Crystal 2:\n  Dataset 1:\n    Name: 1\n"
        f"    Wavelength: {wavelength:.6g}\n"
        "    label #valid %valid min max type\n"
        "    IMEAN 10 100.00% 0 9 J: intensity\n"
        "    SIGIMEAN 10 100.00% 0 9 Q: standard deviation\n"
    )
    executable.write_text(
        "#!/bin/sh\n"
        "cat <<'NASOLVE_HEADER'\n" + header + "NASOLVE_HEADER\n"
        "printf '%s\\n' 'I(+) 1 K' 'SIGI(+) 1 M' 'I(-) 1 K' 'SIGI(-) 1 M'\n"
        "case \"$1\" in\n"
        "  *overall_best_refine_data.mtz) "
        "printf '%s\\n' 'HLAM 1 A' 'HLBM 1 A' 'HLCM 1 A' 'HLDM 1 A';;\n"
        "esac\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_autosol(
    root: Path,
    site_xyz: tuple[float, float, float],
    return_code: int = 0,
    map_files: tuple[str, ...] = (),
) -> Path:
    executable = root / "phenix.autosol"
    ha_record = pdb_record(
        "HETATM",
        1,
        "I",
        "I",
        "A",
        1,
        occupancy=0.93,
        element="I",
        x=site_xyz[0],
        y=site_xyz[1],
        z=site_xyz[2],
    ).rstrip("\n")
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        # Parameter paths from the official AutoSol keyword hierarchy. Keep
        # this external-program contract independent of the command builder.
        "allowed = {'autosol.data', 'autosol.seq_file', 'autosol.labels',\n"
        "           'autosol.atom_type', 'autosol.lambda',\n"
        "           'autosol.phasing.input_partpdb_file',\n"
        "           'autosol.model_building.build',\n"
        "           'autosol.model_building.phase_improve_and_build',\n"
        "           'autosol.general.nproc'}\n"
        "keys = [arg.partition('=')[0] for arg in sys.argv[1:]]\n"
        "if len(keys) != len(allowed) or set(keys) != allowed:\n"
        "    print('Sorry: Ambiguous or unknown AutoSol parameter definition: ' + ', '.join(keys))\n"
        "    raise SystemExit(1)\n"
        "root = Path.cwd() / 'AutoSol_run_1_'\n"
        "root.mkdir()\n"
        "(Path.cwd() / 'received_args.json').write_text(__import__('json').dumps(sys.argv[1:]))\n"
        f"return_code = {return_code}\n"
        "if return_code:\n"
        "    print('AutoSol diagnostic failure')\n"
        "    raise SystemExit(return_code)\n"
        "(root / 'autosol.eff').write_text('sites = None\\nbuild = False\\n"
        "phase_improve_and_build = False\\nnproc = 8\\n')\n"
        f"(root / 'overall_best_ha_pdb.pdb').write_text({ha_record!r} + '\\nEND\\n')\n"
        "(root / 'overall_best_refine_data.mtz').write_bytes(b'mtz')\n"
        f"for name in {map_files!r}:\n"
        "    map_path = root / name\n"
        "    map_path.parent.mkdir(parents=True, exist_ok=True)\n"
        "    map_path.write_bytes(b'density-modified map')\n"
        "print('AutoSol complete')\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_run(root: Path, model_xyz: tuple[float, float, float]) -> Path:
    dataset = root / "dataset"
    frame = root / "MR_frames" / "5W6W"
    run = dataset / "AutoMR" / "run_001"
    phaser = run / "Phaser"
    postmr_model_dir = run / "PostMR" / "Model"
    frame.mkdir(parents=True)
    phaser.mkdir(parents=True)
    postmr_model_dir.mkdir(parents=True)
    input_mtz = dataset / "staraniso_alldata-unique.mtz"
    input_mtz.write_bytes(b"mtz")
    summary = dataset / "summary.html"
    summary.write_text(
        "wavelength [A] = 1.377618\nWavelength               1.37762 A\n"
    )
    frame_model = frame / "C_G.pdb"
    frame_model.write_text("END\n")
    (frame / "seq_base.txt").write_text("GAGC\n\nCTGC\n")
    phaser_model = phaser / "mr_solution.pdb"
    phaser_model.write_text(cryst1() + "END\n")
    prepared_model = postmr_model_dir / "readyset_model.pdb"
    prepared_model.write_text(
        cryst1()
        + pdb_record(
            "HETATM", 1, "I", "C38", "B", 4, element="I",
            x=model_xyz[0], y=model_xyz[1], z=model_xyz[2],
        )
        + "END\n"
    )
    report = {
        "workflow": "automr",
        "stage": "postmr",
        "status": "POSTMR_READY",
        "inputs": {
            "reflections": str(input_mtz.resolve()),
            "reflections_sha256": file_sha256(input_mtz),
            "summary": str(summary.resolve()),
            "model": str(frame_model.resolve()),
        },
        "post_mr_plan": {"sequences": {}},
        "execution": {
            "phaser": {"solution_pdb": str(phaser_model.resolve())},
        },
        "postmr": {
            "prepared_model": str(prepared_model.resolve()),
            "prepared_sha256": file_sha256(prepared_model),
            "anomalous": {
                "autosol_required": True,
                "trigger_elements": ["BR", "I", "SE"],
                "candidates": [{
                    "site": "B:4",
                    "residue": "C38",
                    "atom_name": "I",
                    "alternate": None,
                    "element": "I",
                    "coordinates": list(model_xyz),
                }],
            },
        },
    }
    (run / "report.json").write_text(json.dumps(report))
    return run


class AutoSolTests(unittest.TestCase):
    def test_density_modified_map_reference_survives_relocation(self):
        with symlinked_temporary_directory() as directory:
            root = Path(directory)
            source = root / "source"
            source.mkdir()
            run = make_run(source, (1.0, 2.0, 3.0))
            result = execute_autosol(
                run,
                make_autosol(source, (1.0, 2.0, 3.0), map_files=(
                    "overall_best_denmod_map_coeffs.mtz",
                    "resolve_histograms_1.mtz",
                    "resolve_histograms_2.mtz",
                )),
                make_mtz_dump(source),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_READY")
            outputs = json.loads(result.report_path.read_text())["outputs"]
            reference = outputs["density_modified_map"]
            self.assertEqual(outputs["density_modified_map_status"], "available")
            self.assertEqual(reference["anchor"], "run")
            self.assertEqual(reference["relative_path"],
                "AutoSol/AutoSol_run_1_/overall_best_denmod_map_coeffs.mtz")
            self.assertEqual(reference["sha256"], outputs["density_modified_map_sha256"])
            self.assertEqual(
                json.loads((run / "report.json").read_text())["autosol"]["outputs"],
                outputs,
            )

            relocated = root / "relocated"
            shutil.copytree(source, relocated)
            shutil.rmtree(source)
            relocated_run = relocated / "dataset" / "AutoMR" / "run_001"
            expected = relocated_run / reference["relative_path"]
            self.assertEqual(resolve_artifact_path(reference, relocated_run), expected.resolve())
            self.assertEqual(file_sha256(expected), outputs["density_modified_map_sha256"])
            expected.write_bytes(b"changed map")
            self.assertIsNone(resolve_artifact_path(reference, relocated_run))

    def test_unique_histogram_map_excludes_temporary_products(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            result = execute_autosol(
                run,
                make_autosol(root, (1.0, 2.0, 3.0), map_files=(
                    "resolve_histograms_1.mtz",
                    "PDS/overall_best_denmod_map_coeffs.mtz",
                    "TEMP/overall_best_denmod_map_coeffs.mtz",
                    "TEMP/resolve_histograms_2.mtz",
                )),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_READY")
            outputs = json.loads(result.report_path.read_text())["outputs"]
            self.assertEqual(outputs["density_modified_map_status"], "available")
            self.assertEqual(outputs["density_modified_map"]["relative_path"],
                "AutoSol/AutoSol_run_1_/resolve_histograms_1.mtz")

    def test_missing_density_modified_map_does_not_reject_valid_phases(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            result = execute_autosol(
                run,
                make_autosol(root, (1.0, 2.0, 3.0)),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_READY")
            report = json.loads(result.report_path.read_text())
            self.assertTrue(report["use_for_refinement"])
            outputs = report["outputs"]
            self.assertEqual(outputs["density_modified_map_status"], "missing")
            self.assertIsNone(outputs["density_modified_map"])
            self.assertIsNone(outputs["density_modified_map_sha256"])
            self.assertIn("No density-modified map", outputs["density_modified_map_diagnostic"])
            self.assertTrue(outputs["refinement_data"].endswith("overall_best_refine_data.mtz"))

    def test_ambiguous_density_modified_maps_do_not_reject_valid_phases(self):
        for map_files in (
            (
                "first/overall_best_denmod_map_coeffs.mtz",
                "second/overall_best_denmod_map_coeffs.mtz",
                "resolve_histograms_1.mtz",
            ),
            ("resolve_histograms_1.mtz", "resolve_histograms_2.mtz"),
        ):
            with self.subTest(map_files=map_files), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                run = make_run(root, (1.0, 2.0, 3.0))
                result = execute_autosol(
                    run,
                    make_autosol(root, (1.0, 2.0, 3.0), map_files=map_files),
                    make_mtz_dump(root),
                    environment={"PATH": "/usr/bin:/bin"},
                    processor_count=8,
                )
                self.assertEqual(result.status, "AUTOSOL_READY")
                report = json.loads(result.report_path.read_text())
                self.assertTrue(report["use_for_refinement"])
                outputs = report["outputs"]
                self.assertEqual(outputs["density_modified_map_status"], "ambiguous")
                self.assertIsNone(outputs["density_modified_map"])
                self.assertIsNone(outputs["density_modified_map_sha256"])
                self.assertIn("found 2", outputs["density_modified_map_diagnostic"])

    def test_unreadable_density_modified_map_does_not_reject_valid_phases(self):
        def checksum(path):
            if path.name == "overall_best_denmod_map_coeffs.mtz":
                raise OSError("map unavailable")
            return file_sha256(path)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            with patch("nasolve.autosol.file_sha256", side_effect=checksum):
                result = execute_autosol(
                    run,
                    make_autosol(root, (1.0, 2.0, 3.0), map_files=(
                        "overall_best_denmod_map_coeffs.mtz",
                    )),
                    make_mtz_dump(root),
                    environment={"PATH": "/usr/bin:/bin"},
                    processor_count=8,
                )
            self.assertEqual(result.status, "AUTOSOL_READY")
            outputs = json.loads(result.report_path.read_text())["outputs"]
            self.assertEqual(outputs["density_modified_map_status"], "unavailable")
            self.assertIsNone(outputs["density_modified_map"])
            self.assertIsNone(outputs["density_modified_map_sha256"])
            self.assertIn("map unavailable", outputs["density_modified_map_diagnostic"])

    def test_anchored_frame_sequence_requires_and_matches_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            report = json.loads((run / "report.json").read_text())
            sequence = run / "Model" / "seq_base.txt"
            sequence.parent.mkdir(exist_ok=True)
            sequence.write_text("GAGC\n\nCTGC\n")
            reference = artifact_reference(sequence, run)
            report["inputs"]["frame_sequence"] = reference

            with self.assertRaisesRegex(
                AutoSolPreparationError,
                "no checksum",
            ):
                _sequence_source(report, run)

            reference["sha256"] = 123
            with self.assertRaisesRegex(
                AutoSolPreparationError,
                "malformed",
            ):
                _sequence_source(report, run)

            reference["sha256"] = file_sha256(sequence)
            selected, _ = _sequence_source(report, run)
            self.assertEqual(selected, sequence.resolve())
            sequence.write_text("changed\n")
            with self.assertRaisesRegex(
                AutoSolPreparationError,
                "does not match",
            ):
                _sequence_source(report, run)

    def test_windows_model_path_does_not_become_bare_sequence_lookup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            report = json.loads((run / "report.json").read_text())
            report["inputs"]["model"] = (
                r"C:\source\MR_frames\5W6W\E_G.pdb"
            )
            report["post_mr_plan"]["sequences"] = {
                "A": "GAGC",
                "B": "CTGC",
            }
            (run / "seq_base.txt").write_text("UNRELATED\n")

            selected, generated = _sequence_source(report, run)

            self.assertIsNone(selected)
            self.assertEqual(generated, "GAGC\n\nCTGC\n")

    def test_legacy_report_without_frozen_checksums_still_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            del report["inputs"]["reflections_sha256"]
            del report["postmr"]["prepared_sha256"]
            report_path.write_text(json.dumps(report))

            result = execute_autosol(
                run,
                make_autosol(root, (1.0, 2.0, 3.0)),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )

            result_report = json.loads(result.report_path.read_text())
            self.assertEqual(
                result.status,
                "AUTOSOL_READY",
                result_report.get("failure_reason"),
            )
            self.assertTrue(
                (result.autosol_directory / "received_args.json").is_file()
            )

    def test_legacy_frame_sequence_rebases_without_exact_pair_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            checkout = root / "checkout"
            run = make_run(source, (1.0, 2.0, 3.0))
            (source / "pyproject.toml").write_text("[project]\nname='test'\n")
            report = json.loads((run / "report.json").read_text())
            missing_pair_model = source / "MR_frames" / "5W6W" / "E_G.pdb"
            missing_pair_model.write_text("END\n")
            report["inputs"]["model"] = str(missing_pair_model.resolve())
            report["post_mr_plan"]["sequences"] = {}
            report["frame"] = {"catalogue_directory": "5W6W"}
            (run / "report.json").write_text(json.dumps(report))
            shutil.copytree(source, checkout)
            (checkout / "MR_frames" / "5W6W" / "E_G.pdb").unlink()
            relocated_run = checkout / "dataset" / "AutoMR" / "run_001"
            relocated_report = json.loads(
                (relocated_run / "report.json").read_text()
            )

            sequence, description = _sequence_source(
                relocated_report,
                relocated_run,
            )

            self.assertEqual(
                sequence,
                (checkout / "MR_frames" / "5W6W" / "seq_base.txt").resolve(),
            )
            self.assertIn("frame", str(description))
            shutil.rmtree(source)
            sequence_without_source, _ = _sequence_source(
                relocated_report,
                relocated_run,
            )
            self.assertEqual(
                sequence_without_source,
                (checkout / "MR_frames" / "5W6W" / "seq_base.txt").resolve(),
            )

    def test_changed_frozen_reflections_stop_before_autosol_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            report = json.loads((run / "report.json").read_text())
            Path(report["inputs"]["reflections"]).write_bytes(b"changed")

            result = execute_autosol(
                run,
                make_autosol(root, (1.0, 2.0, 3.0)),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
            )

            self.assertEqual(result.status, "AUTOSOL_WARNING")
            warning = json.loads(result.report_path.read_text())
            self.assertIn("reflections", warning["failure_reason"])
            self.assertFalse((result.autosol_directory / "received_args.json").exists())

    def test_changed_prepared_model_stops_before_autosol_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            report = json.loads((run / "report.json").read_text())
            Path(report["postmr"]["prepared_model"]).write_text("changed\n")

            result = execute_autosol(
                run,
                make_autosol(root, (1.0, 2.0, 3.0)),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
            )

            self.assertEqual(result.status, "AUTOSOL_WARNING")
            warning = json.loads(result.report_path.read_text())
            self.assertIn("prepared model", warning["failure_reason"])
            self.assertFalse((result.autosol_directory / "received_args.json").exists())

    def test_malformed_present_checksum_stops_before_autosol_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root, (1.0, 2.0, 3.0))
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["inputs"]["reflections_sha256"] = 123
            report_path.write_text(json.dumps(report))

            result = execute_autosol(
                run,
                make_autosol(root, (1.0, 2.0, 3.0)),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
            )

            self.assertEqual(result.status, "AUTOSOL_WARNING")
            warning = json.loads(result.report_path.read_text())
            self.assertIn("malformed", warning["failure_reason"])
            self.assertFalse((result.autosol_directory / "received_args.json").exists())

    def test_wavelength_prefers_precise_assignment(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = Path(directory) / "summary.html"
            summary.write_text(
                "wavelength [A] = 1.377618\nWavelength 1.37762 A\n"
            )
            value, source = read_wavelength(summary)
            self.assertEqual(value, 1.377618)
            self.assertIn("high-precision", source)

    def test_conflicting_wavelengths_stop(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = Path(directory) / "summary.html"
            summary.write_text("wavelength [A] = 1.377618\nWavelength 0.91990 A\n")
            with self.assertRaisesRegex(AutoSolPreparationError, "Conflicting"):
                read_wavelength(summary)

    def test_r3_symmetry_match_finds_third_operator(self):
        model = hex_cart(0.2, 0.1, 0.3)
        transformed = hex_cart(-0.1, -0.2, 0.3)
        site = (transformed[0] + 1.0, transformed[1], transformed[2])
        result = nearest_symmetry_distance(
            model, site, (50.0, 50.0, 60.0, 90.0, 90.0, 120.0), "R 3"
        )
        self.assertAlmostEqual(result["distance"], 1.0, places=6)
        self.assertEqual(result["operator"], "-x+y,-x,z")

    def test_execute_is_mr_guided_never_builds_and_validates_site(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = hex_cart(0.2, 0.1, 0.3)
            transformed = hex_cart(-0.1, -0.2, 0.3)
            site = (transformed[0] + 1.0, transformed[1], transformed[2])
            run = make_run(root, model)
            result = execute_autosol(
                run,
                make_autosol(root, site),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_READY")
            self.assertEqual(result.matched_distance, 1.0)
            args = json.loads((result.autosol_directory / "received_args.json").read_text())
            self.assertIn("autosol.model_building.build=False", args)
            self.assertIn("autosol.model_building.phase_improve_and_build=False", args)
            self.assertIn("autosol.general.nproc=8", args)
            self.assertFalse(any(argument.partition("=")[0].split(".")[-1] == "sites" for argument in args))
            model_argument = next(arg for arg in args if arg.startswith("autosol.phasing.input_partpdb_file="))
            self.assertTrue(model_argument.endswith("/Phaser/mr_solution.pdb"))
            report = json.loads((run / "report.json").read_text())
            self.assertEqual(report["stage"], "autosol")
            self.assertFalse(report["autosol"]["effective_parameters"]["build"])
            outputs = report["autosol"]["outputs"]
            self.assertRegex(outputs["heavy_atom_model_sha256"], r"^[0-9a-f]{64}$")
            self.assertRegex(outputs["refinement_data_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(
                report["autosol"]["site_validation"]["accepted"][0]
                ["nearest_site"]["operator"],
                "-x+y,-x,z",
            )

    def test_fully_qualified_parameters_preserve_mr_sad_inputs_and_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            run = make_run(root, (1.0, 2.0, 3.0))
            result = execute_autosol(
                run,
                make_autosol(root, (1.0, 2.0, 3.0)),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_READY", result.message)
            args = json.loads((result.autosol_directory / "received_args.json").read_text())
            self.assertEqual(dict(arg.split("=", 1) for arg in args), {
                "autosol.data": str(root / "dataset/staraniso_alldata-unique.mtz"),
                "autosol.seq_file": str(run / "AutoSol/sequence_input.txt"),
                "autosol.labels": "I(+) SIGI(+) I(-) SIGI(-)",
                "autosol.atom_type": "I",
                "autosol.lambda": "1.37762",
                "autosol.phasing.input_partpdb_file": str(run / "Phaser/mr_solution.pdb"),
                "autosol.model_building.build": "False",
                "autosol.model_building.phase_improve_and_build": "False",
                "autosol.general.nproc": "8",
            })
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["command"][1:], args)

    def test_distant_ha_sites_stop_before_report_acceptance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = hex_cart(0.2, 0.1, 0.3)
            run = make_run(root, model)
            result = execute_autosol(
                run,
                make_autosol(root, (13.0, 14.0, 15.0), map_files=(
                    "overall_best_denmod_map_coeffs.mtz",
                )),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_WARNING")
            report = json.loads((run / "report.json").read_text())
            self.assertFalse(report["autosol"]["use_for_refinement"])
            self.assertIn("continue without AutoSol phases", report["autosol"]["continuation"])
            self.assertEqual(
                report["autosol"]["outputs"]["density_modified_map_status"], "available"
            )

    def test_intermediate_distance_requests_review_but_continues(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = hex_cart(0.2, 0.1, 0.3)
            transformed = hex_cart(-0.1, -0.2, 0.3)
            site = (transformed[0] + 6.0, transformed[1], transformed[2])
            run = make_run(root, model)
            result = execute_autosol(
                run,
                make_autosol(root, site, map_files=("resolve_histograms_1.mtz",)),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_REVIEW")
            self.assertAlmostEqual(result.matched_distance, 6.0)
            report = json.loads((run / "report.json").read_text())
            self.assertFalse(report["autosol"]["use_for_refinement"])
            self.assertEqual(
                len(report["autosol"]["site_validation"]["review_candidates"]), 1
            )
            self.assertEqual(
                report["autosol"]["outputs"]["density_modified_map_status"], "available"
            )

    def test_autosol_program_failure_warns_and_continues(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = hex_cart(0.2, 0.1, 0.3)
            run = make_run(root, model)
            result = execute_autosol(
                run,
                make_autosol(root, (0.0, 0.0, 0.0), return_code=5),
                make_mtz_dump(root),
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_WARNING")
            self.assertIsNone(result.refinement_data)
            report = json.loads((run / "report.json").read_text())
            self.assertEqual(report["autosol"]["status"], "AUTOSOL_WARNING")
            self.assertIn("status 5", report["autosol"]["failure_reason"])


if __name__ == "__main__":
    unittest.main()
