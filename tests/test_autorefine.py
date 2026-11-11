import json
import shutil
import stat
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import nasolve.checkpoints as checkpoint_module
from nasolve.autorefine import (
    AutoRefineError,
    DATA_MANAGER_FILE_SCOPED,
    audit_terminal_phosphate_geometry,
    LEGACY_EXPLICIT,
    _cached_file_referencer,
    _run_report,
    anomalous_selections,
    with_postmr_metal_sites,
    build_reflection_plan,
    execute_autorefine,
    parse_refinement_statistics,
    reflection_selector_policy,
    write_terminal_phosphate_protection,
)
from nasolve.checkpoints import initialize_registry, list_checkpoints
from nasolve.model_assessment import file_sha256

from .helpers import pdb_record


PHENIX_120 = "1.20.1-4487"
PHENIX_21 = "2.1-6048"


def make_refine_run(root: Path, *, autosol: bool = True) -> Path:
    dataset = root / "dataset"
    run = dataset / "AutoMR" / "run_001"
    model_dir = run / "PostMR" / "Model"
    restraints_dir = run / "PostMR" / "Restraints"
    readyset_dir = run / "PostMR" / "ReadySet"
    model_dir.mkdir(parents=True)
    restraints_dir.mkdir(parents=True)
    readyset_dir.mkdir(parents=True)
    observations = dataset / "staraniso.mtz"
    observations.parent.mkdir(exist_ok=True)
    observations.write_bytes(b"observations")
    model = model_dir / "readyset_model.pdb"
    model.write_text(
        pdb_record("ATOM", 1, "P", "S6G", "A", 12)
        + pdb_record("HETATM", 2, "I", "C38", "B", 4, element="I")
        + "END\n"
    )
    phil = restraints_dir / "narestraints.phil"
    phil.write_text("geometry_restraints {}\n")
    readyset_cif = readyset_dir / "prepared_model.ligands.cif"
    readyset_cif.write_text("data_components\n")
    report = {
        "workflow": "automr",
        "stage": "postmr",
        "status": "POSTMR_READY",
        "inputs": {"reflections": str(observations.resolve())},
        "postmr": {
            "status": "POSTMR_READY",
            "created_utc": "2026-08-28T12:00:00+00:00",
            "prepared_model": str(model.resolve()),
            "restraints": [str(phil.resolve())],
            "readyset": {"generated_ligand_cif": str(readyset_cif.resolve())},
            "mutation_actions": [
                {"site": "A:12", "before": "DC", "after": "S6G"},
                {"site": "B:4", "before": "DG", "after": "C38"},
            ],
            "anomalous": {
                "candidates": [{
                    "site": "B:4",
                    "residue": "C38",
                    "atom_name": "I",
                    "alternate": None,
                    "element": "I",
                    "coordinates": [1.0, 2.0, 3.0],
                }],
            },
        },
    }
    if autosol:
        phase = run / "AutoSol" / "overall_best_refine_data.mtz"
        phase.parent.mkdir(parents=True)
        phase.write_bytes(b"phases")
        report["stage"] = "autosol"
        report["status"] = "AUTOSOL_READY"
        report["autosol"] = {
            "status": "AUTOSOL_READY",
            "use_for_refinement": True,
            "wavelength": 1.377618,
            "outputs": {
                "refinement_data": str(phase.resolve()),
                "refinement_data_sha256": file_sha256(phase),
            },
        }
    (run / "report.json").write_text(json.dumps(report))
    return run


def make_mtz_dump(root: Path, *, anomalous: bool = True) -> Path:
    executable = root / "phenix.mtz.dump"
    anomalous_lines = (
        "printf '%s\\n' 'F(+) 1 G' 'SIGF(+) 1 L' 'F(-) 1 G' 'SIGF(-) 1 L'\n"
        if anomalous else ""
    )
    executable.write_text(
        "#!/bin/sh\n"
        "case \"$1\" in\n"
        "  *overall_best*) printf '%s\\n' 'HLAM 1 A' 'HLBM 1 A' 'HLCM 1 A' 'HLDM 1 A';;\n"
        "  *) printf '%s\\n' 'IMEAN 1 J' 'SIGIMEAN 1 Q' 'FreeR_flag 1 I';\n"
        f"     {anomalous_lines}"
        "     printf '%s\\n' 'Resolution range: 20.0 5.27';\n"
        "     ;;\n"
        "esac\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable



def make_pdb_interpretation(root: Path) -> Path:
    executable = root / "phenix.pdb_interpretation"

    def bond(first: str, second: str, ideal: float, sigma: float) -> str:
        return (
            f'bond pdb=" {first} DA  D   1 "\n'
            f'     pdb=" {second} DA  D   1 "\n'
            "  ideal  model  delta    sigma   weight residual\n"
            f"  {ideal:.3f}  {ideal:.3f} 0.000 {sigma:.2e} 1.00e+00 0.00e+00\n"
        )

    def angle(
        first: str,
        second: str,
        third: str,
        ideal: float,
        sigma: float,
    ) -> str:
        return (
            f'angle pdb=" {first} DA  D   1 "\n'
            f'      pdb=" {second} DA  D   1 "\n'
            f'      pdb=" {third} DA  D   1 "\n'
            "  ideal  model  delta    sigma   weight residual\n"
            f"  {ideal:.2f}  {ideal:.2f} 0.00 {sigma:.2e} 1.00e+00 0.00e+00\n"
        )

    geometry = "".join([
        bond("P", "O5'", 1.593, 0.010),
        bond("P", "OP1", 1.480, 0.020),
        bond("P", "OP2", 1.480, 0.020),
        bond("P", "OP3", 1.480, 0.020),
        angle("OP1", "P", "OP2", 120.00, 3.00),
        angle("OP1", "P", "O5'", 109.00, 3.00),
        angle("OP2", "P", "O5'", 108.00, 3.00),
        angle("OP1", "P", "OP3", 109.47, 3.00),
        angle("OP2", "P", "OP3", 109.47, 3.00),
        angle("O5'", "P", "OP3", 109.47, 3.00),
        angle("P", "O5'", "C5'", 120.90, 1.60),
    ])
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        "if 'write_geo=True' not in sys.argv:\n"
        "    print('ERROR: expected write_geo=True')\n"
        "    raise SystemExit(9)\n"
        "model = next(Path(arg) for arg in sys.argv[1:] if arg.endswith('.pdb'))\n"
        f"geometry = {geometry!r}\n"
        "Path(model.name + '.geo').write_text(geometry)\n"
        "print('Wrote pre-refinement geometry')\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_refine(
    root: Path,
    *,
    final_work: float = 0.274,
    final_free: float = 0.301,
    return_code: int = 0,
    preflight_return_code: int = 0,
    legacy_single_mtz: bool = False,
    phase_change: str | None = None,
    reject_data_manager: bool = False,
) -> Path:
    executable = root / "phenix.refine"
    mtz_outputs = (
        "Path('refined_001.mtz').write_bytes(b'legacy maps and reflections')\n"
        if legacy_single_mtz
        else (
            "Path('refined_001_map_coeffs.mtz').write_bytes(b'maps')\n"
            "Path('refined_001_data.mtz').write_bytes(b'refine evidence')\n"
        )
    )
    phase_change_code = {
        None: "",
        "modify": (
            "phase = next(Path(arg) for arg in sys.argv[1:] "
            "if arg.endswith('overall_best_refine_data.mtz'))\n"
            "phase.write_bytes(b'changed while refining')\n"
        ),
        "delete": (
            "phase = next(Path(arg) for arg in sys.argv[1:] "
            "if arg.endswith('overall_best_refine_data.mtz'))\n"
            "phase.unlink()\n"
        ),
    }.get(phase_change)
    if phase_change_code is None:
        raise ValueError(f"Unsupported phase test change: {phase_change}")
    script = (
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import json, shutil, sys\n"
        "if '--dry_run' in sys.argv:\n"
        f"    preflight_return_code = {preflight_return_code}\n"
        f"    reject_data_manager = {reject_data_manager!r}\n"
        "    parameter_files = [Path(arg) for arg in sys.argv[1:] if arg.endswith('.params')]\n"
        "    if reject_data_manager and any('data_manager {' in path.read_text() for path in parameter_files):\n"
        "        print('ERROR: Unused parameter definitions: data_manager.miller_array.file')\n"
        "        raise SystemExit(8)\n"
        "    print('Dry-run parameters accepted' if not preflight_return_code else 'ERROR bad parameter')\n"
        "    raise SystemExit(preflight_return_code)\n"
        "Path('received_args.json').write_text(json.dumps(sys.argv[1:]))\n"
        f"return_code = {return_code}\n"
        "if return_code:\n"
        "    print('Sorry: deliberate refinement failure')\n"
        "    raise SystemExit(return_code)\n"
        "model = next(Path(arg) for arg in sys.argv[1:] if arg.endswith('.pdb'))\n"
        "shutil.copyfile(model, 'refined_001.pdb')\n"
        "Path('refined_001.cif').write_text('data_model\\n')\n"
        "Path('refined_001.reflections.cif').write_text('data_reflections\\n')\n"
    ) + mtz_outputs + phase_change_code + (
        "print('Start: r_work = 0.412 r_free = 0.438 bonds = 0.012 angles = 1.82')\n"
        "print('Macrocycle 1: r_work = 0.331 r_free = 0.359')\n"
        f"print('Final: r_work = {final_work} r_free = {final_free} bonds = 0.009 angles = 1.21')\n"
        "print('clashscore = 4.8')\n"
        "print('Anomalous scatterer group:')\n"
        "print('  Selection: \\\"chain B and resid 4 and name I\\\"')\n"
        "print('  Number of selected scatterers: 1')\n"
        "print('  f_prime: -1.8')\n"
        "print('  f_double_prime: 7.91907')\n"
    )
    executable.write_text(script)
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    phenix_python = root / "phenix.python"
    phenix_python.write_text(
        "#!/bin/sh\n"
        "case \"$*\" in\n"
        "  *audit_free_r.py*) printf '%s\\n' 'NASOLVE_FREE_R_AUDIT_JSON:{\"array_anomalous\": true, \"stored_observations\": 2017, \"independent_friedel_groups\": 1010, \"paired_friedel_groups\": 1007, \"free_independent_groups\": 45, \"free_fraction\": 0.044554, \"inconsistent_friedel_flag_groups\": 0, \"test_flag_value\": 0, \"resolution_shells\": [{\"groups\": 101, \"free_groups\": 4}, {\"groups\": 101, \"free_groups\": 5}]}';;\n"
        "  *NASOLVE_SCATTERING_JSON*) printf '%s\\n' 'NASOLVE_SCATTERING_JSON:{\"I\": [-1.8, 5.75]}';;\n"
        "  *) printf '%s\\n' '5.750252723693848';;\n"
        "esac\n"
    )
    phenix_python.chmod(phenix_python.stat().st_mode | stat.S_IXUSR)
    return executable


class AutoRefineTests(unittest.TestCase):
    def test_reflection_selector_policy_is_explicit_and_fail_closed(self):
        self.assertEqual(
            reflection_selector_policy(PHENIX_120).mode,
            LEGACY_EXPLICIT,
        )
        self.assertEqual(
            reflection_selector_policy(PHENIX_21).mode,
            DATA_MANAGER_FILE_SCOPED,
        )
        for version in (None, "unknown", "2.x", "2.0-0000", "2.3-9999", "3.0-0000"):
            with self.subTest(version=version):
                with self.assertRaises(AutoRefineError):
                    reflection_selector_policy(version)  # type: ignore[arg-type]

    def test_phenix_22_uses_file_scoped_selector_policy(self):
        for version in ("2.2", "2.2.1", "2.2.1-6174"):
            with self.subTest(version=version):
                self.assertEqual(
                    reflection_selector_policy(version).mode,
                    DATA_MANAGER_FILE_SCOPED,
                )

    def test_phenix_22_scopes_observations_flags_and_phases_to_their_files(self):
        for autosol in (False, True):
            with self.subTest(autosol=autosol), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                run = make_refine_run(root, autosol=autosol)
                report_path = run / "report.json"
                report = json.loads(report_path.read_text())
                if not autosol:
                    report["postmr"]["anomalous"]["candidates"] = []
                    report_path.write_text(json.dumps(report))

                result = execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version="2.2.1-6174",
                    environment={"PATH": "/usr/bin:/bin"},
                )

                self.assertEqual(result.status, "AUTOREFINE_READY")
                payload = json.loads(result.report_path.read_text())
                self.assertEqual(payload["phenix_version"], "2.2.1-6174")
                self.assertEqual(
                    payload["reflection_selector_mode"], DATA_MANAGER_FILE_SCOPED,
                )
                params = (result.round_directory / "autorefine.params").read_text()
                arrays = params.split("  miller_array {")[1:]
                self.assertEqual(len(arrays), 2 if autosol else 1)
                observations = (root / "dataset" / "staraniso.mtz").resolve()
                self.assertIn(f'file = "{observations}"', arrays[0])
                labels = "F(+),SIGF(+),F(-),SIGF(-)" if autosol else "IMEAN,SIGIMEAN"
                self.assertIn(f'name = "{labels}"', arrays[0])
                self.assertIn('name = "FreeR_flag"', arrays[0])
                self.assertNotIn('name = "HLAM,HLBM,HLCM,HLDM"', arrays[0])
                if autosol:
                    phases = (run / "AutoSol" / "overall_best_refine_data.mtz").resolve()
                    self.assertIn(f'file = "{phases}"', arrays[1])
                    self.assertIn('name = "HLAM,HLBM,HLCM,HLDM"', arrays[1])
                    self.assertNotIn('name = "FreeR_flag"', arrays[1])
                    self.assertNotIn(f'name = "{labels}"', arrays[1])

    def test_geometry_conservative_recipe_reduces_xray_weight_scale(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root),
                make_mtz_dump(root),
                phenix_version="2.2.1-6174",
                environment={"PATH": "/usr/bin:/bin"},
                recipe="AutoRefine/geometry-conservative",
                macro_cycles=1,
            )

            params = (result.round_directory / "autorefine.params").read_text()
            self.assertIn("optimize_xyz_weight = True", params)
            self.assertIn("wxc_scale = 0.1", params)
            self.assertIn("write_final_geo_file = True", params)

            payload = json.loads(result.report_path.read_text())
            self.assertEqual(payload["recipe"], "AutoRefine/geometry-conservative")

    def test_terminal_phosphate_geometry_audit_uses_phenix_native_sigmas(self):
        with tempfile.TemporaryDirectory() as directory:
            geo = Path(directory) / "final.geo"

            def bond(a, b, ideal=1.480, model=1.480, sigma=0.020):
                return (
                    f'bond pdb="{a:>4} DC  D   1 "\n'
                    f'     pdb="{b:>4} DC  D   1 "\n'
                    "  ideal  model  delta    sigma   weight residual\n"
                    f"  {ideal:.3f}  {model:.3f} {ideal-model:6.3f} {sigma:.2e} 1.00e+00 0.00e+00\n"
                )

            def angle(a, b, c, ideal, model=None, sigma=3.0):
                model = ideal if model is None else model
                return (
                    f'angle pdb="{a:>4} DC  D   1 "\n'
                    f'      pdb="{b:>4} DC  D   1 "\n'
                    f'      pdb="{c:>4} DC  D   1 "\n'
                    "    ideal   model   delta    sigma   weight residual\n"
                    f"   {ideal:.2f}  {model:.2f}  {ideal-model:.2f} {sigma:.2e} 1.00e+00 0.00e+00\n"
                )

            clean = "".join([
                bond("P", "O5'", 1.593, 1.593, 0.010),
                bond("P", "OP1"),
                bond("P", "OP2"),
                bond("P", "OP3"),
                angle("OP1", "P", "OP2", 120.00),
                angle("OP1", "P", "OP3", 109.47),
                angle("OP2", "P", "OP3", 109.47),
                angle("O5'", "P", "OP1", 109.00),
                angle("O5'", "P", "OP2", 108.00),
                angle("O5'", "P", "OP3", 109.47),
                angle("P", "O5'", "C5'", 120.90, sigma=1.60),
            ])
            geo.write_text(clean)

            audit = audit_terminal_phosphate_geometry(geo, ("D:1",))
            self.assertEqual(audit["status"], "PASS")
            self.assertFalse(audit["requires_review"])
            self.assertEqual(audit["sites"][0]["restraint_count"], 11)

            geo.write_text(
                clean.replace(
                    angle("OP1", "P", "OP3", 109.47),
                    angle("OP1", "P", "OP3", 109.47, model=88.47),
                )
            )
            audit = audit_terminal_phosphate_geometry(geo, ("D:1",))
            self.assertEqual(audit["status"], "REVIEW")
            self.assertTrue(audit["requires_review"])
            self.assertAlmostEqual(
                audit["sites"][0]["max_sigma_deviation"],
                7.0,
                places=6,
            )

    def test_terminal_phosphate_protection_can_target_clean_declared_site(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "terminal_phosphate_protection.phil"
            audit = {
                "status": "PASS",
                "requires_review": False,
                "sites": [{
                    "site": "D:1",
                    "requires_review": False,
                    "restraints": [
                        {"kind": "angle", "atoms": ["OP1", "P", "OP2"], "ideal": 120.00},
                        {"kind": "angle", "atoms": ["OP1", "P", "OP3"], "ideal": 109.47},
                        {"kind": "angle", "atoms": ["OP2", "P", "OP3"], "ideal": 109.47},
                        {"kind": "angle", "atoms": ["O5'", "P", "OP1"], "ideal": 109.00},
                        {"kind": "angle", "atoms": ["O5'", "P", "OP2"], "ideal": 108.00},
                        {"kind": "angle", "atoms": ["O5'", "P", "OP3"], "ideal": 109.47},
                    ],
                }],
            }

            record = write_terminal_phosphate_protection(
                audit,
                destination,
                protect_sites=("D:1",),
            )

            self.assertEqual(record["schema_version"], 1)
            self.assertEqual(record["kind"], "terminal-phosphate-protection")
            self.assertEqual(record["sites"], ["D:1"])
            self.assertEqual(record["angle_count"], 6)
            self.assertEqual(
                record["ideal_source"],
                "source-checkpoint-final-phenix-geometry-audit",
            )
            self.assertEqual(destination.read_text().count("action = *change"), 6)

            with self.assertRaises(
                AutoRefineError,
                msg="PASS audit must not trigger Doctor-style protection implicitly",
            ):
                write_terminal_phosphate_protection(
                    audit,
                    Path(directory) / "implicit.phil",
                )

    def test_terminal_phosphate_protection_fails_closed_for_missing_declared_site(self):
        with tempfile.TemporaryDirectory() as directory:
            audit = {
                "status": "PASS",
                "requires_review": False,
                "sites": [{
                    "site": "D:1",
                    "requires_review": False,
                    "restraints": [
                        {"kind": "angle", "atoms": ["OP1", "P", "OP2"], "ideal": 120.00},
                        {"kind": "angle", "atoms": ["OP1", "P", "OP3"], "ideal": 109.47},
                        {"kind": "angle", "atoms": ["OP2", "P", "OP3"], "ideal": 109.47},
                        {"kind": "angle", "atoms": ["O5'", "P", "OP1"], "ideal": 109.00},
                        {"kind": "angle", "atoms": ["O5'", "P", "OP2"], "ideal": 108.00},
                        {"kind": "angle", "atoms": ["O5'", "P", "OP3"], "ideal": 109.47},
                    ],
                }],
            }

            with self.assertRaises(
                AutoRefineError,
                msg="Explicit proactive protection must fail if the requested site is absent",
            ):
                write_terminal_phosphate_protection(
                    audit,
                    Path(directory) / "missing.phil",
                    protect_sites=("A:1",),
                )

    def test_terminal_geometry_review_blocks_numerical_auto_promotion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            forced_audit = {
                "status": "REVIEW",
                "requires_review": True,
                "sigma_threshold": 5.0,
                "sites": [{"site": "D:1", "max_sigma_deviation": 7.0}],
            }
            with patch(
                "nasolve.autorefine.audit_terminal_phosphate_geometry",
                return_value=forced_audit,
            ):
                result = execute_autorefine(
                    run,
                    make_refine(root, final_work=0.20, final_free=0.24),
                    make_mtz_dump(root),
                    phenix_version="2.2.1-6174",
                    environment={"PATH": "/usr/bin:/bin"},
                    macro_cycles=1,
                )

            self.assertEqual(result.status, "AUTOREFINE_REVIEW")
            self.assertFalse(result.selected_as_current)
            payload = json.loads(result.report_path.read_text())
            self.assertTrue(
                payload["acceptance"]["terminal_geometry_review"]
            )

    def test_extra_restraint_is_passed_to_refine_and_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            extra = root / "terminal-protection.phil"
            extra.write_text("refinement.geometry_restraints.edits {}\n")

            result = execute_autorefine(
                run,
                make_refine(root, final_work=0.244, final_free=0.267),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
                macro_cycles=1,
                extra_restraints=(extra,),
                auto_select_success=False,
            )

            payload = json.loads(result.report_path.read_text())
            expected = str(extra.resolve())
            self.assertEqual(payload["inputs"]["extra_restraints"], [expected])
            self.assertIn(expected, payload["inputs"]["restraints"])
            self.assertIn(expected, payload["command"])

    def test_terminal_protection_metadata_is_frozen_and_inherited(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            extra = run / "terminal-protection.phil"
            extra.write_text("refinement.geometry_restraints.edits {}\n")
            protection = {
                "schema_version": 1,
                "kind": "terminal-phosphate-protection",
                "path": str(extra.resolve()),
                "sites": ["D:1"],
                "angle_count": 6,
                "sigma": 1.0,
                "ideal_source": "test-phenix-geometry",
                "mechanism": "phenix-action-change",
            }

            first = execute_autorefine(
                run,
                make_refine(root, final_work=0.244, final_free=0.267),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
                macro_cycles=1,
                extra_restraints=(extra,),
                terminal_geometry_protection=protection,
                auto_select_success=False,
            )

            registry = json.loads(
                (run / "AutoRefine" / "checkpoints.json").read_text()
            )
            first_checkpoint = next(
                item for item in registry["checkpoints"]
                if item["id"] == first.checkpoint_id
            )
            frozen = first_checkpoint["terminal_geometry_protection"]
            self.assertNotIn("path", frozen)
            self.assertEqual(frozen["sites"], ["D:1"])
            self.assertEqual(frozen["angle_count"], 6)
            self.assertEqual(frozen["ideal_source"], "test-phenix-geometry")
            self.assertEqual(frozen["restraint"]["anchor"], "run")
            self.assertEqual(frozen["restraint"]["sha256"], file_sha256(extra))

            second = execute_autorefine(
                run,
                make_refine(root, final_work=0.243, final_free=0.266),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
                from_checkpoint=first.checkpoint_id,
                macro_cycles=1,
                auto_select_success=False,
            )
            registry = json.loads(
                (run / "AutoRefine" / "checkpoints.json").read_text()
            )
            second_checkpoint = next(
                item for item in registry["checkpoints"]
                if item["id"] == second.checkpoint_id
            )
            self.assertEqual(
                second_checkpoint["terminal_geometry_protection"],
                frozen,
            )

            with self.assertRaisesRegex(
                AutoRefineError,
                "already carries terminal-geometry protection",
            ):
                execute_autorefine(
                    run,
                    make_refine(root, final_work=0.242, final_free=0.265),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                    from_checkpoint=first.checkpoint_id,
                    macro_cycles=1,
                    extra_restraints=(extra,),
                    terminal_geometry_protection=protection,
                    auto_select_success=False,
                )


    def test_declared_terminal_phosphate_is_protected_before_first_refinement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            make_pdb_interpretation(root)
            refine = make_refine(root, final_work=0.244, final_free=0.267)
            mtz_dump = make_mtz_dump(root)

            with (
                patch(
                    "nasolve.autorefine.requested_op3_sites",
                    return_value=("D:1",),
                ),
                patch(
                    "nasolve.autorefine.validate_refined_model",
                    return_value={"validated": True},
                ),
            ):
                first = execute_autorefine(
                    run,
                    refine,
                    mtz_dump,
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                    macro_cycles=1,
                    auto_select_success=False,
                )

            protection = (
                first.round_directory / "terminal_phosphate_protection.phil"
            )
            source_geometry = (
                first.round_directory / "terminal_geometry_pre_refine.geo"
            )
            source_log = (
                first.round_directory / "terminal_geometry_pre_refine.log"
            )
            self.assertTrue(protection.is_file())
            self.assertTrue(source_geometry.is_file())
            self.assertTrue(source_log.is_file())
            self.assertEqual(
                protection.read_text().count("action = *change"),
                6,
            )

            payload = json.loads(first.report_path.read_text())
            expected_protection = str(protection.resolve())
            self.assertIn(
                expected_protection,
                payload["inputs"]["restraints"],
            )
            self.assertIn(
                expected_protection,
                payload["command"],
            )

            registry = json.loads(
                (run / "AutoRefine" / "checkpoints.json").read_text()
            )
            first_checkpoint = next(
                item for item in registry["checkpoints"]
                if item["id"] == first.checkpoint_id
            )
            semantic = first_checkpoint["terminal_geometry_protection"]
            self.assertEqual(semantic["sites"], ["D:1"])
            self.assertEqual(semantic["angle_count"], 6)
            self.assertEqual(
                semantic["ideal_source"],
                "pre-refinement-phenix.pdb_interpretation-geometry",
            )
            self.assertEqual(semantic["restraint"]["anchor"], "run")
            self.assertEqual(semantic["source_geometry"]["anchor"], "run")
            self.assertEqual(semantic["source_log"]["anchor"], "run")

            # An inherited protected child must not need to invoke interpretation
            # again and must carry exactly one protection artifact.
            (root / "phenix.pdb_interpretation").unlink()
            with (
                patch(
                    "nasolve.autorefine.requested_op3_sites",
                    return_value=("D:1",),
                ),
                patch(
                    "nasolve.autorefine.validate_refined_model",
                    return_value={"validated": True},
                ),
            ):
                second = execute_autorefine(
                    run,
                    make_refine(root, final_work=0.243, final_free=0.266),
                    mtz_dump,
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                    from_checkpoint=first.checkpoint_id,
                    macro_cycles=1,
                    auto_select_success=False,
                )

            second_payload = json.loads(second.report_path.read_text())
            protection_inputs = [
                value for value in second_payload["inputs"]["restraints"]
                if value.endswith("terminal_phosphate_protection.phil")
            ]
            self.assertEqual(protection_inputs, [expected_protection])

            registry = json.loads(
                (run / "AutoRefine" / "checkpoints.json").read_text()
            )
            second_checkpoint = next(
                item for item in registry["checkpoints"]
                if item["id"] == second.checkpoint_id
            )
            self.assertEqual(
                second_checkpoint["terminal_geometry_protection"],
                semantic,
            )

    def test_missing_interpreter_fails_before_numbered_round_is_created(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)

            with (
                patch(
                    "nasolve.autorefine.requested_op3_sites",
                    return_value=("D:1",),
                ),
                patch(
                    "nasolve.autorefine.validate_refined_model",
                    return_value={"validated": True},
                ),
                self.assertRaisesRegex(
                    AutoRefineError,
                    "phenix.pdb_interpretation was not found",
                ),
            ):
                execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                    macro_cycles=1,
                    auto_select_success=False,
                )

            self.assertFalse(
                (run / "AutoRefine" / "round_001").exists()
            )

    def test_phenix_22_preflight_failure_prevents_refinement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root, preflight_return_code=1),
                make_mtz_dump(root),
                phenix_version="2.2.1-6174",
                environment={"PATH": "/usr/bin:/bin"},
            )
            self.assertEqual(result.status, "AUTOREFINE_FAILED")
            self.assertFalse((result.round_directory / "received_args.json").exists())

    def test_unknown_phenix_version_fails_before_creating_refinement_state(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)

            with self.assertRaisesRegex(AutoRefineError, "will not guess"):
                execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version="unknown",
                    environment={"PATH": "/usr/bin:/bin"},
                )

            self.assertFalse((run / "AutoRefine").exists())

    def test_non_object_report_is_rejected_cleanly(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            (run / "report.json").write_text(json.dumps([]))

            with self.assertRaisesRegex(AutoRefineError, "not a JSON object"):
                _run_report(run)

    def test_phenix_120_keeps_exact_anomalous_phase_selectors_without_data_manager(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root, reject_data_manager=True),
                make_mtz_dump(root),
                phenix_version=PHENIX_120,
                environment={"PATH": "/usr/bin:/bin"},
            )

            self.assertEqual(result.status, "AUTOREFINE_READY")
            payload = json.loads(result.report_path.read_text())
            params = (result.round_directory / "autorefine.params").read_text()
            args = json.loads((result.round_directory / "received_args.json").read_text())
            run_report = json.loads((run / "report.json").read_text())
            observations = run_report["inputs"]["reflections"]
            phases = run_report["autosol"]["outputs"]["refinement_data"]
            expected = {
                f"xray_data.file_name={observations}",
                "xray_data.labels=F(+),SIGF(+),F(-),SIGF(-)",
                f"xray_data.r_free_flags.file_name={observations}",
                "xray_data.r_free_flags.label=FreeR_flag",
                "xray_data.r_free_flags.test_flag_value=0",
                "xray_data.r_free_flags.generate=False",
                f"experimental_phases.file_name={phases}",
                "experimental_phases.labels=HLAM,HLBM,HLCM,HLDM",
                'main.target="auto"',
                "xray_data.force_anomalous_flag_to_be_equal_to=True",
            }
            self.assertTrue(expected.issubset(args))
            self.assertNotIn("data_manager", params)
            self.assertNotIn("miller_array", params)
            self.assertNotIn("--unused_ok", payload["command"])
            self.assertNotIn("--unused_ok", payload["preflight_command"])
            self.assertEqual(payload["phenix_version"], PHENIX_120)
            self.assertEqual(payload["reflection_selector_mode"], LEGACY_EXPLICIT)
            self.assertEqual(run_report["autorefine"]["phenix_version"], PHENIX_120)
            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            checkpoint = registry["checkpoints"][-1]
            self.assertEqual(checkpoint["phenix_version"], PHENIX_120)
            self.assertEqual(checkpoint["reflection_selector_mode"], LEGACY_EXPLICIT)

    def test_phenix_120_keeps_exact_mean_selectors_without_data_manager(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["postmr"]["anomalous"]["candidates"] = []
            report_path.write_text(json.dumps(report))

            result = execute_autorefine(
                run,
                make_refine(root, reject_data_manager=True),
                make_mtz_dump(root),
                phenix_version=PHENIX_120,
                environment={"PATH": "/usr/bin:/bin"},
            )

            params = (result.round_directory / "autorefine.params").read_text()
            args = json.loads((result.round_directory / "received_args.json").read_text())
            observations = json.loads((run / "report.json").read_text())["inputs"]["reflections"]
            self.assertNotIn("data_manager", params)
            self.assertIn(f"xray_data.file_name={observations}", args)
            self.assertIn("xray_data.labels=IMEAN,SIGIMEAN", args)
            self.assertIn(f"xray_data.r_free_flags.file_name={observations}", args)
            self.assertIn("xray_data.r_free_flags.label=FreeR_flag", args)
            self.assertIn("xray_data.r_free_flags.test_flag_value=0", args)
            self.assertIn("xray_data.r_free_flags.generate=False", args)
            self.assertIn("xray_data.force_anomalous_flag_to_be_equal_to=False", args)
            self.assertFalse(any(arg.startswith("experimental_phases.") for arg in args))

    def test_plan_uses_friedel_amplitudes_phases_and_exact_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            report = json.loads((run / "report.json").read_text())
            observations = Path(report["inputs"]["reflections"])
            phase = Path(report["autosol"]["outputs"]["refinement_data"])
            plan = build_reflection_plan(
                report, observations, make_mtz_dump(root), phase_file=phase
            )
            self.assertEqual(plan.observation_labels, ("F(+)", "SIGF(+)", "F(-)", "SIGF(-)"))
            self.assertEqual(plan.phase_labels, ("HLAM", "HLBM", "HLCM", "HLDM"))
            self.assertTrue(plan.anomalous)
            self.assertEqual(
                anomalous_selections(report),
                ("chain B and resid 4 and name I",),
            )

    def test_postmr_metal_sites_are_refined_as_anomalous_scatterers(self):
        # Regression: only nucleotide I/Br/Se candidates reached phenix.refine.
        report = {"postmr": {"anomalous": {
            "candidates": [{"site": "A:8", "atom_name": "I5", "element": "I"}],
            "metal_sites": [{"site": "Z:1", "atom_name": "AG", "element": "AG", "kind": "metal"}],
        }}}
        self.assertEqual(
            anomalous_selections(report),
            ("chain A and resid 8 and name I5", "chain Z and resid 1 and name AG"),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            metal_only = json.loads((run / "report.json").read_text())
            metal_only["postmr"]["anomalous"] = {
                "candidates": [], "autosol_required": False,
                "metal_sites": [{"site": "Z:1", "atom_name": "AG", "element": "AG"}],
            }
            self.assertEqual(anomalous_selections(metal_only), ("chain Z and resid 1 and name AG",))
            observations = Path(metal_only["inputs"]["reflections"])
            plan = build_reflection_plan(metal_only, observations, make_mtz_dump(root))
            self.assertTrue(plan.anomalous)
            self.assertEqual(plan.observation_labels, ("F(+)", "SIGF(+)", "F(-)", "SIGF(-)"))

    def test_postmr_from_before_metal_detection_gets_metal_sites_from_its_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            report = json.loads((run / "report.json").read_text())
            report["postmr"]["anomalous"].pop("metal_sites", None)
            prepared = run / "prepared_with_silver.pdb"
            prepared.write_text(
                "HETATM    1 AG    AG B 101       1.000   2.000   3.000  1.00 20.00          AG\nEND\n"
            )
            report["postmr"]["prepared_model"] = str(prepared)
            report["postmr"]["prepared_sha256"] = file_sha256(prepared)
            derived = with_postmr_metal_sites(report, run)
            self.assertEqual(
                [(s["site"], s["atom_name"]) for s in derived["postmr"]["anomalous"]["metal_sites"]],
                [("B:101", "AG")],
            )
            self.assertIn("chain B and resid 101 and name AG", anomalous_selections(derived))
            self.assertNotIn("metal_sites", report["postmr"]["anomalous"])  # not rewritten
            # A newer PostMR record (key present, even empty) is used as is.
            report["postmr"]["anomalous"]["metal_sites"] = []
            self.assertIs(with_postmr_metal_sites(report, run), report)
            prepared.write_text("changed\n")
            del report["postmr"]["anomalous"]["metal_sites"]
            with self.assertRaisesRegex(AutoRefineError, "changed after PostMR"):
                with_postmr_metal_sites(report, run)

    def test_success_is_quiet_checkpointed_and_reused_by_next_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            refine = make_refine(root)
            mtz_dump = make_mtz_dump(root)
            progress: list[tuple[str, Path]] = []
            first = execute_autorefine(
                run,
                refine,
                mtz_dump,
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
                progress=lambda checkpoint, log: progress.append((checkpoint, log)),
            )
            self.assertEqual(first.status, "AUTOREFINE_READY")
            self.assertEqual(progress[0][0], "refine-001")
            self.assertEqual(first.statistics["r_work"], 0.274)
            self.assertEqual(first.statistics["clashscore"], 4.8)
            anomalous = first.statistics["anomalous_scatterers"][0]
            self.assertAlmostEqual(anomalous["refined_f_double_prime"], 7.91907)
            self.assertAlmostEqual(anomalous["calculated_f_double_prime"], 5.750252723693848)
            self.assertAlmostEqual(anomalous["model_occupancy"], 1.0)
            self.assertAlmostEqual(anomalous["apparent_anomalous_occupancy"], 1.377168, places=5)
            self.assertAlmostEqual(anomalous["resolution_limit"], 5.27)
            self.assertTrue(first.log_path.is_file())
            self.assertEqual(first.model_cif.name, "refined_001.cif")
            self.assertEqual(first.reflection_cif.name, "refined_001.reflections.cif")
            self.assertTrue((first.round_directory / "metrics.tsv").is_file())
            args = json.loads((first.round_directory / "received_args.json").read_text())
            self.assertIn('main.target="auto"', args)
            self.assertIn("xray_data.labels=F(+),SIGF(+),F(-),SIGF(-)", args)
            self.assertFalse(
                any(
                    argument.startswith("data_manager.miller_array.labels.name=")
                    for argument in args
                )
            )
            run_report = json.loads((run / "report.json").read_text())
            self.assertIn(
                "xray_data.file_name=" + run_report["inputs"]["reflections"],
                args,
            )
            params = (first.round_directory / "autorefine.params").read_text()
            self.assertEqual(params.count("  miller_array {"), 2)
            self.assertIn(
                f'file = "{run_report["inputs"]["reflections"]}"',
                params,
            )
            self.assertIn(
                f'file = "{run_report["autosol"]["outputs"]["refinement_data"]}"',
                params,
            )
            self.assertIn('name = "F(+),SIGF(+),F(-),SIGF(-)"', params)
            self.assertIn('name = "FreeR_flag"', params)
            self.assertIn('name = "HLAM,HLBM,HLCM,HLDM"', params)
            self.assertIn('target = "auto"', params)
            self.assertIn("strategy = *individual_sites *individual_sites_real_space", params)
            self.assertIn("individual = all", params)
            self.assertIn("group_adp_refinement_mode = one_adp_group_per_residue", params)
            self.assertIn("group = all", params)
            self.assertIn('selection = "chain B and resid 4 and name I"', params)
            self.assertIn("group_anomalous", params)
            observation_file_marker = f'file = "{run_report["inputs"]["reflections"]}"'
            phase_file_marker = (
                f'file = "{run_report["autosol"]["outputs"]["refinement_data"]}"'
            )
            observation_label_marker = 'name = "F(+),SIGF(+),F(-),SIGF(-)"'
            self.assertLess(
                params.index(observation_file_marker),
                params.index(observation_label_marker),
            )
            self.assertLess(
                params.index(observation_label_marker),
                params.index('name = "FreeR_flag"'),
            )
            self.assertLess(params.index('name = "FreeR_flag"'), params.index(phase_file_marker))
            self.assertLess(
                params.index(phase_file_marker),
                params.index('name = "HLAM,HLBM,HLCM,HLDM"'),
            )
            first_report = json.loads(first.report_path.read_text())
            self.assertEqual(first_report["phenix_version"], PHENIX_21)
            self.assertEqual(
                first_report["reflection_selector_mode"],
                DATA_MANAGER_FILE_SCOPED,
            )
            first_registry = json.loads(
                (run / "AutoRefine" / "checkpoints.json").read_text()
            )
            first_checkpoint = first_registry["checkpoints"][-1]
            self.assertEqual(first_checkpoint["phenix_version"], PHENIX_21)
            self.assertEqual(
                first_checkpoint["reflection_selector_mode"],
                DATA_MANAGER_FILE_SCOPED,
            )

            second = execute_autorefine(
                run,
                refine,
                mtz_dump,
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            self.assertEqual(second.checkpoint_id, "refine-002")
            self.assertEqual(second.parent_checkpoint, "refine-001")
            second_report = json.loads(second.report_path.read_text())
            self.assertEqual(second_report["inputs"]["model"], str(first.model_path))
            self.assertTrue(second_report["inputs"]["authoritative_observations"].endswith("staraniso.mtz"))

    def test_review_checkpoint_does_not_replace_current(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root, final_work=0.31, final_free=0.34),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
            )
            self.assertEqual(result.status, "AUTOREFINE_REVIEW")
            records, current, _ = list_checkpoints(run)
            self.assertEqual(current, "postmr")
            self.assertEqual(records[-1].status, "REVIEW")
            self.assertTrue(records[-1].usable)

    def test_controlled_trial_can_disable_phases_real_space_occupancy_and_anomalous(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
                use_experimental_phases=False,
                real_space_sites=False,
                adp_mode="group",
                refine_occupancies=False,
                anomalous_mode="off",
                auto_select_success=False,
            )
            payload = json.loads(result.report_path.read_text())
            params = (result.round_directory / "autorefine.params").read_text()
            self.assertEqual(payload["inputs"]["phase_labels"], [])
            self.assertNotIn("individual_sites_real_space", params)
            self.assertNotIn("occupancies", params)
            self.assertNotIn("group_anomalous", params)
            self.assertFalse(result.selected_as_current)

    def test_missing_friedel_amplitudes_after_autosol_is_final_error_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root),
                make_mtz_dump(root, anomalous=False),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
            )
            self.assertEqual(result.status, "AUTOREFINE_ANOMALOUS_FALLBACK")
            self.assertEqual(result.exit_code, 2)
            payload = json.loads(result.report_path.read_text())
            self.assertFalse(payload["refinement"]["anomalous"])
            self.assertEqual(payload["inputs"]["observation_labels"], ["IMEAN", "SIGIMEAN"])
            self.assertNotIn("group_anomalous", payload["refinement"]["strategies"])

    def test_mean_observations_use_data_manager_selector_without_autosol(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["postmr"]["anomalous"]["candidates"] = []
            report_path.write_text(json.dumps(report))

            result = execute_autorefine(
                run,
                make_refine(root),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
            )

            self.assertEqual(result.status, "AUTOREFINE_READY")
            args = json.loads((result.round_directory / "received_args.json").read_text())
            self.assertIn("xray_data.labels=IMEAN,SIGIMEAN", args)
            self.assertFalse(
                any(
                    argument.startswith("data_manager.miller_array.labels.name=")
                    for argument in args
                )
            )
            params = (result.round_directory / "autorefine.params").read_text()
            self.assertEqual(params.count("  miller_array {"), 1)
            self.assertIn('name = "IMEAN,SIGIMEAN"', params)
            self.assertIn('name = "FreeR_flag"', params)
            self.assertNotIn("HLAM,HLBM,HLCM,HLDM", params)

    def test_existing_pre_autosol_registry_recovers_phases_after_move(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_run = make_refine_run(root / "collaborator", autosol=False)
            initialize_registry(source_run)
            phase = source_run / "AutoSol" / "overall_best_refine_data.mtz"
            phase.parent.mkdir(parents=True)
            phase.write_bytes(b"phases")
            report_path = source_run / "report.json"
            report = json.loads(report_path.read_text())
            report["autosol"] = {
                "status": "AUTOSOL_READY",
                "use_for_refinement": True,
                "outputs": {
                    "refinement_data": str(phase.resolve()),
                    "refinement_data_sha256": file_sha256(phase),
                },
            }
            report_path.write_text(json.dumps(report))
            checkout_dataset = root / "checkout" / "dataset"
            shutil.copytree(source_run.parent.parent, checkout_dataset)
            shutil.rmtree(root / "collaborator")
            run = checkout_dataset / "AutoMR" / source_run.name

            with patch(
                "nasolve.autorefine.file_sha256", wraps=file_sha256
            ) as hash_file:
                result = execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                )

            params = (result.round_directory / "autorefine.params").read_text()
            self.assertIn(str(run.resolve() / "AutoSol"), params)
            self.assertIn("HLAM,HLBM,HLCM,HLDM", params)
            phase_path = (run / "AutoSol" / phase.name).resolve()
            phase_hashes = [
                call
                for call in hash_file.call_args_list
                if Path(call.args[0]).resolve() == phase_path
            ]
            self.assertEqual(len(phase_hashes), 1)

    def test_inherited_phase_is_hashed_once_and_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            phase = (run / "AutoSol" / "overall_best_refine_data.mtz").resolve()
            initialize_registry(run)

            with (
                patch(
                    "nasolve.autorefine.file_sha256", wraps=file_sha256
                ) as hash_file,
                patch(
                    "nasolve.checkpoints._sha256",
                    wraps=checkpoint_module._sha256,
                ) as checkpoint_hash,
            ):
                result = execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                )

            phase_hashes = [
                call
                for call in hash_file.call_args_list
                if Path(call.args[0]).resolve() == phase
            ]
            checkpoint_phase_hashes = [
                call
                for call in checkpoint_hash.call_args_list
                if Path(call.args[0]).resolve() == phase
            ]
            self.assertEqual(
                len(phase_hashes) + len(checkpoint_phase_hashes),
                1,
            )
            registry = json.loads(
                (run / "AutoRefine" / "checkpoints.json").read_text()
            )
            parent, child = registry["checkpoints"]
            self.assertEqual(child["id"], result.checkpoint_id)
            self.assertEqual(child["phases"]["sha256"], parent["phases"]["sha256"])
            self.assertEqual(
                child["phases"]["labels"],
                ["HLAM", "HLBM", "HLCM", "HLDM"],
            )

    def test_inherited_phase_checksum_drift_still_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            initialize_registry(run)
            phase = run / "AutoSol" / "overall_best_refine_data.mtz"
            phase.write_bytes(b"changed phases")

            with self.assertRaisesRegex(
                AutoRefineError, "changed after it was frozen"
            ):
                execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                )

    def test_phase_changed_during_phenix_creates_only_a_failed_checkpoint(self):
        for phase_change in ("modify", "delete"):
            with self.subTest(phase_change=phase_change):
                with tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    run = make_refine_run(root)
                    initialize_registry(run)

                    result = execute_autorefine(
                        run,
                        make_refine(root, phase_change=phase_change),
                        make_mtz_dump(root),
                        phenix_version=PHENIX_21,
                        environment={"PATH": "/usr/bin:/bin"},
                    )

                    registry = json.loads(
                        (run / "AutoRefine" / "checkpoints.json").read_text()
                    )
                    self.assertEqual(result.status, "AUTOREFINE_FAILED")
                    self.assertIn("while Phenix was running", result.message)
                    self.assertEqual(registry["current"], "postmr")
                    self.assertEqual(len(registry["checkpoints"]), 2)
                    self.assertFalse(registry["checkpoints"][-1]["usable"])

                    phase = run / "AutoSol" / "overall_best_refine_data.mtz"
                    phase.write_bytes(b"phases")
                    retry = execute_autorefine(
                        run,
                        make_refine(root),
                        make_mtz_dump(root),
                        phenix_version=PHENIX_21,
                        environment={"PATH": "/usr/bin:/bin"},
                    )
                    self.assertEqual(retry.checkpoint_id, "refine-002")
                    self.assertEqual(retry.status, "AUTOREFINE_READY")

    def test_checkpoint_outputs_hash_each_unique_file_once(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)

            with patch(
                "nasolve.autorefine.file_sha256", wraps=file_sha256
            ) as hash_file:
                result = execute_autorefine(
                    run,
                    make_refine(root, legacy_single_mtz=True),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                )

            primary_mtz = (result.round_directory / "refined_001.mtz").resolve()
            primary_hashes = [
                call
                for call in hash_file.call_args_list
                if Path(call.args[0]).resolve() == primary_mtz
            ]
            self.assertEqual(len(primary_hashes), 1)
            registry = json.loads(
                (run / "AutoRefine" / "checkpoints.json").read_text()
            )
            outputs = registry["checkpoints"][-1]["outputs"]
            self.assertEqual(
                outputs["map_coefficients"], outputs["refinement_reflections"]
            )

    def test_cached_output_references_are_independent_copies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root / "dataset" / "AutoMR" / "run_001"
            artifact = run / "AutoRefine" / "shared.mtz"
            artifact.parent.mkdir(parents=True)
            artifact.write_bytes(b"shared output")

            with patch(
                "nasolve.autorefine.file_sha256", wraps=file_sha256
            ) as hash_file:
                reference = _cached_file_referencer(run)
                first = reference(artifact)
                first["consumer"] = "maps"
                second = reference(artifact)

            self.assertNotIn("consumer", second)
            self.assertEqual(len(hash_file.call_args_list), 1)

    def test_approved_missing_phase_file_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            initialize_registry(run)
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["autosol"] = {
                "status": "AUTOSOL_READY",
                "use_for_refinement": True,
                "outputs": {
                    "refinement_data": str(run / "AutoSol" / "missing.mtz")
                },
            }
            report_path.write_text(json.dumps(report))

            with self.assertRaisesRegex(AutoRefineError, "phase data is missing"):
                execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                )

    def test_approved_phase_checksum_drift_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            initialize_registry(run)
            phase = run / "AutoSol" / "overall_best_refine_data.mtz"
            phase.parent.mkdir(parents=True)
            phase.write_bytes(b"original phases")
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["autosol"] = {
                "status": "AUTOSOL_READY",
                "use_for_refinement": True,
                "outputs": {
                    "refinement_data": str(phase),
                    "refinement_data_sha256": file_sha256(phase),
                },
            }
            report_path.write_text(json.dumps(report))
            phase.write_bytes(b"changed phases")

            with self.assertRaisesRegex(AutoRefineError, "checksum validation"):
                execute_autorefine(
                    run,
                    make_refine(root),
                    make_mtz_dump(root),
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                )

    def test_failed_phenix_round_is_visible_and_not_reusable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root, return_code=7),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
            )
            self.assertEqual(result.status, "AUTOREFINE_FAILED")
            records, current, _ = list_checkpoints(run)
            self.assertEqual(current, "postmr")
            self.assertFalse(records[-1].usable)
            self.assertEqual(records[-1].status, "FAILED")

    def test_failed_preflight_never_starts_real_refinement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = execute_autorefine(
                run,
                make_refine(root, preflight_return_code=9),
                make_mtz_dump(root),
                phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"},
            )
            self.assertEqual(result.status, "AUTOREFINE_FAILED")
            self.assertFalse((result.round_directory / "received_args.json").exists())
            self.assertIn("preflight failed", result.log_path.read_text().lower())
            payload = json.loads(result.report_path.read_text())
            self.assertEqual(payload["phenix_version"], PHENIX_21)
            self.assertEqual(payload["reflection_selector_mode"], DATA_MANAGER_FILE_SCOPED)
            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            failed = registry["checkpoints"][-1]
            self.assertEqual(failed["status"], "FAILED")
            self.assertEqual(failed["phenix_version"], PHENIX_21)
            self.assertEqual(failed["reflection_selector_mode"], DATA_MANAGER_FILE_SCOPED)

    def test_statistics_parser_keeps_cycle_history(self):
        stats = parse_refinement_statistics(
            "Start: r_work=0.40 r_free=0.44 bonds=0.01 angles=1.5\n"
            "cycle: r_work=0.35 r_free=0.39\n"
            "Final: r_work=0.28 r_free=0.31 bonds=0.009 angles=1.2\n"
            "clashscore: 5.1\n"
            "WARNING geometry needs inspection\n"
        )
        self.assertEqual(len(stats["cycle_series"]), 3)
        self.assertEqual(stats["r_free_minus_r_work"], 0.03)
        self.assertEqual(stats["bond_rmsd"], 0.009)
        self.assertEqual(stats["diagnostics"], ["WARNING geometry needs inspection"])

    def test_missing_free_r_is_rejected_without_generation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            dump = root / "phenix.mtz.dump"
            dump.write_text("#!/bin/sh\necho 'IMEAN 1 J'\necho 'SIGIMEAN 1 Q'\n")
            dump.chmod(dump.stat().st_mode | stat.S_IXUSR)
            with self.assertRaisesRegex(AutoRefineError, "Free-R"):
                execute_autorefine(
                    run,
                    make_refine(root),
                    dump,
                    phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"},
                )


if __name__ == "__main__":
    unittest.main()
