"""AutoFDP from a model and an MTZ: preflight, step 1 (B start) and step 2
(one f''-only round). No NASolve run, report or checkpoint is involved. Step 3
(occupancy adjustment) is in test_autofdp_occupancy.py; here the stand-in f''
does not depend on occupancy.

phenix.refine, phenix.mtz.dump and phenix.python are simulated. The
stand-in phenix.refine reads the anomalous groups NASolve wrote into its
parameters and echoes them, with a chosen refined f'', in phenix.refine's log
format, so the round, the f'' comparison and the record run end to end.
"""

import contextlib
import io
import json
import shlex
import stat
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from nasolve.autofdp import (
    STATUS_CONVERGED,
    STATUS_FAILED,
    STATUS_REVIEW,
    STATUS_ROUND_DONE,
    AutoFDPError,
    execute_autofdp,
    held_parameter_changes,
    parse_theory_overrides,
    plan_autofdp,
    read_reflections,
)
from nasolve.cli import main
from nasolve.config import AppConfig
from nasolve.model_assessment import file_sha256
from nasolve.phenix_runtime import PhenixInstallation

from .helpers import mtz_dump_columns_text, pdb_record
from .test_autofdp_bstart import bridged_model_text
from .test_autorefine import PHENIX_21, make_mtz_dump

ENV = {"PATH": "/usr/bin:/bin"}
SILVER = "chain B and resid 101 and name AG"
MERCURY = "chain C and resid 1 and name HG"


def make_inputs(root: Path, *, mercury: bool = False) -> tuple[Path, Path]:
    """A refined model (Ag B:101, occupancy 0.90, B 379.99, bridging A:12 DC
    with B 60 and B:4 DU with B 40) and an MTZ in a dataset folder."""
    text = bridged_model_text()
    if mercury:  # a second metal beside A:20 DG (B 5), away from the Ag
        text = text.replace("END\n", pdb_record(
            "HETATM", 15, "HG", "HG", "C", 1, element="HG", x=22.0, y=21.0, z=20.0,
            b_factor=150.0, occupancy=0.5,
        ) + "END\n")
    dataset = root / "7CU_Ag"
    round_dir = dataset / "AutoMR" / "run_001" / "AutoRefine" / "round_002"
    round_dir.mkdir(parents=True)
    model = round_dir / "refined_001.pdb"
    model.write_text(text)
    mtz = dataset / "staraniso_alldata-unique.mtz"
    mtz.write_bytes(b"reflections")
    return model, mtz


def make_fdp_refine(root: Path, fdp: dict[str, float] | None = None, *,
                    change_b: bool = False) -> Path:
    """Stand-in phenix.refine that echoes each anomalous group it was given,
    with f'' from ``fdp`` (groups not in ``fdp`` are not reported)."""
    fdp = {SILVER: 2.6, MERCURY: 9.0} if fdp is None else fdp
    executable = root / "phenix.refine"
    edit = (
        "lines = Path('refined_001.pdb').read_text().splitlines(keepends=True)\n"
        "lines[0] = lines[0][:60] + '%6.2f' % (float(lines[0][60:66]) + 5) + lines[0][66:]\n"
        "Path('refined_001.pdb').write_text(''.join(lines))\n"
        if change_b else ""
    )
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, re, shutil, sys\n"
        "from pathlib import Path\n"
        "if '--dry_run' in sys.argv:\n"
        "    print('Dry-run parameters accepted')\n"
        "    raise SystemExit(0)\n"
        "Path('received_args.json').write_text(json.dumps(sys.argv[1:]))\n"
        "params = next(Path(a) for a in sys.argv[1:] if a.endswith('.params')).read_text()\n"
        "model = next(Path(a) for a in sys.argv[1:] if a.endswith('.pdb'))\n"
        "shutil.copyfile(model, 'refined_001.pdb')\n"
        + edit +
        "Path('refined_001.cif').write_text('data_model\\n')\n"
        "Path('refined_001_map_coeffs.mtz').write_bytes(b'maps')\n"
        "Path('refined_001_data.mtz').write_bytes(b'refine evidence')\n"
        "print('Start: r_work = 0.2165 r_free = 0.2124 bonds = 0.012 angles = 1.82')\n"
        "print('Final: r_work = 0.2161 r_free = 0.2121 bonds = 0.012 angles = 1.82')\n"
        f"fdp = {fdp!r}\n"
        "for block in re.findall(r'group \\{(.*?)\\}', params, re.S):\n"
        "    selection = re.search(r'selection = \"(.*?)\"', block).group(1)\n"
        "    if selection not in fdp:\n"
        "        continue\n"
        "    print('Anomalous scatterer group:')\n"
        "    print(f'  Selection: \"{selection}\"')\n"
        "    print('  Number of selected scatterers: 1')\n"
        "    print('  f_prime: ' + re.search(r'f_prime = (\\S+)', block).group(1))\n"
        "    print(f'  f_double_prime: {fdp[selection]}')\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    python = root / "phenix.python"
    python.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' 'NASOLVE_SCATTERING_JSON:{\"AG\": [-0.12, 2.03], \"HG\": [-12.5, 10.1]}'\n"
    )
    python.chmod(python.stat().st_mode | stat.S_IXUSR)
    return executable


def make_ready_set(root: Path, *, mode: str = "copy", generated: tuple[str, ...] = ()) -> Path:
    """Stand-in phenix.ready_set. ``copy`` writes <stem>.updated.pdb unchanged
    (and <stem>.ligands.cif for the ``generated`` residue codes); ``noop`` says
    "No unknown residues" and writes nothing; ``fail``, ``change_b`` and
    ``add_h`` misbehave."""
    executable = root / "phenix.ready_set"
    cif = ""
    if generated:
        cif = "data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.group\n" + "".join(
            f"{code} non-polymer\n" for code in generated
        ) + "".join(
            f"\ndata_comp_{code}\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n"
            f"_chem_comp_atom.type_symbol\n{code} C1 C\n{code} C2 C\n"
            for code in generated
        )
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, sys\n"
        "from pathlib import Path\n"
        "Path('readyset_args.json').write_text(json.dumps(sys.argv[1:]))\n"
        f"mode = {mode!r}\n"
        "model = Path(sys.argv[1])\n"
        "if mode == 'fail':\n"
        "    print('Sorry: cannot process'); raise SystemExit(1)\n"
        "if mode == 'noop':\n"
        "    print('No unknown residues found'); raise SystemExit(0)\n"
        "lines = model.read_text().splitlines(keepends=True)\n"
        "if mode == 'change_b':\n"
        "    lines[0] = lines[0][:60] + '%6.2f' % (float(lines[0][60:66]) + 9) + lines[0][66:]\n"
        "if mode == 'add_h':\n"
        "    lines.insert(0, lines[0][:12] + ' H1 ' + lines[0][16:76] + ' H\\n')\n"
        "Path(model.stem + '.updated.pdb').write_text(''.join(lines))\n"
        f"cif = {cif!r}\n"
        "if cif:\n"
        "    Path(model.stem + '.ligands.cif').write_text(cif)\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def custom_mtz_dump(root: Path, text: str) -> Path:
    executable = root / "custom.mtz.dump"
    executable.write_text("#!/bin/sh\ncat <<'EOF'\n" + text + "EOF\n")
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def _plan(root, model, mtz, **keywords):
    return plan_autofdp(
        model, mtz, keywords.pop("refine", None) or make_fdp_refine(root),
        keywords.pop("dump", None) or make_mtz_dump(root, wavelength=1.00743),
        environment=ENV, **keywords,
    )


def _execute(plan, root, **keywords):
    return execute_autofdp(
        plan, keywords.pop("refine", None) or make_fdp_refine(root),
        keywords.pop("ready_set", None) or make_ready_set(root),
        phenix_version=PHENIX_21, environment=ENV, processor_count=1, **keywords,
    )


class PreflightTests(unittest.TestCase):
    def test_everything_is_read_from_the_model_and_the_mtz(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            plan = _plan(root, model, mtz)
            reflections = plan.reflections
            self.assertEqual(reflections.observation_labels, ("F(+)", "SIGF(+)", "F(-)", "SIGF(-)"))
            self.assertEqual((reflections.free_r_label, reflections.free_r_test_value), ("FreeR_flag", 0))
            self.assertEqual(reflections.wavelength, 1.00743)
            self.assertEqual(reflections.wavelength_source["method"], "MTZ header (phenix.mtz.dump)")
            self.assertEqual(plan.theory["AG"]["f_double_prime"], 2.03)
            self.assertEqual(plan.theory["AG"]["source"], "cctbx-henke")
            self.assertEqual(plan.selections, (SILVER,))
            self.assertEqual(plan.b_start.sites[0].b_after, 50.0)
            self.assertEqual(plan.output_root, (mtz.parent / "AutoFDP").resolve())
            self.assertFalse((mtz.parent / "AutoFDP").exists())  # the plan writes nothing

    def test_mtz_requirements(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            merged_only = custom_mtz_dump(root, mtz_dump_columns_text(
                (("H", "H"), ("K", "H"), ("L", "H"), ("F", "F"), ("SIGF", "Q"),
                 ("FreeR_flag", "I")), wavelength=1.00743,
            ))
            with self.assertRaisesRegex(AutoFDPError, "no anomalous data"):
                _plan(root, model, mtz, dump=merged_only)
            no_wavelength = make_mtz_dump(root)
            with self.assertRaisesRegex(AutoFDPError, "give the wavelength with --wavelength"):
                _plan(root, model, mtz, dump=no_wavelength)
            given = read_reflections(mtz, no_wavelength, ENV, wavelength=1.0074)
            self.assertEqual((given.wavelength, given.wavelength_source["method"]), (1.0074, "--wavelength"))
            differs = read_reflections(mtz, make_mtz_dump(root, wavelength=1.00743), ENV, wavelength=0.9795)
            self.assertIn("differs from the MTZ header", differs.wavelength_source["warning"])
            with self.assertRaisesRegex(AutoFDPError, "between 0.1 and 5"):
                read_reflections(mtz, no_wavelength, ENV, wavelength=12.0)

            phenix_flags = custom_mtz_dump(root, mtz_dump_columns_text(
                (("H", "H"), ("K", "H"), ("L", "H"), ("F(+)", "G"), ("SIGF(+)", "L"),
                 ("F(-)", "G"), ("SIGF(-)", "L"), ("R-free-flags", "I")), wavelength=1.00743,
            ))
            self.assertEqual(read_reflections(mtz, phenix_flags, ENV).free_r_test_value, 1)
            no_flags = custom_mtz_dump(root, mtz_dump_columns_text(
                (("H", "H"), ("K", "H"), ("L", "H"), ("F(+)", "G"), ("SIGF(+)", "L"),
                 ("F(-)", "G"), ("SIGF(-)", "L")), wavelength=1.00743,
            ))
            with self.assertRaisesRegex(AutoFDPError, "no recognized Free-R label"):
                read_reflections(mtz, no_flags, ENV)
            with self.assertRaisesRegex(AutoFDPError, "No MTZ file"):
                read_reflections(root / "missing.mtz", no_flags, ENV)

    def test_theory_override_parsing(self):
        self.assertEqual(parse_theory_overrides(["hg=10.17", "AG=2.06"]), {"HG": 10.17, "AG": 2.06})
        for bad in (["HG"], ["HG=x"], ["HG=-1"], ["H1=2"], ["HG=1", "HG=2"]):
            with self.subTest(bad=bad), self.assertRaises(AutoFDPError):
                parse_theory_overrides(bad)


class FdpRoundTests(unittest.TestCase):
    def test_round_refines_only_the_anomalous_scatterers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            original = model.read_bytes()
            result = _execute(_plan(root, model, mtz), root)
            # Within tolerance in round_001: converged without adjustment.
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)
            self.assertEqual(result.rounds, (result.round,))
            self.assertEqual(result.round.status, STATUS_ROUND_DONE)
            self.assertEqual(result.adjustments, ())
            run = mtz.parent.resolve() / "AutoFDP" / "run_001"
            self.assertEqual(result.run_directory, run)
            self.assertEqual(model.read_bytes(), original)   # the given model is never changed
            self.assertEqual((run / "source_model.pdb").read_bytes(), original)

            # Every strategy off except group_anomalous; f' fixed, f'' from 0.
            params = (run / "round_001" / "autofdp.params").read_text()
            self.assertIn("strategy = *group_anomalous\n", params)
            self.assertIn(f'selection = "{SILVER}"', params)
            self.assertIn("f_prime = -0.12", params)
            self.assertIn("f_double_prime = 0.0", params)
            self.assertIn("refine = f_prime *f_double_prime", params)
            self.assertIn("optimize_xyz_weight = False", params)
            self.assertIn("optimize_adp_weight = False", params)
            self.assertIn("number_of_macro_cycles = 3", params)
            received = json.loads((run / "round_001" / "received_args.json").read_text())
            self.assertIn(str(run / "readyset_model.pdb"), received)
            # ReadySet ran on the step-1 model without hydrogens; the fixture's
            # DU is a curated residue, so NASolve's DU dictionary went with it.
            readyset_args = json.loads((run / "ReadySet" / "readyset_args.json").read_text())
            self.assertEqual(readyset_args[0], str(run / "start_model.pdb"))
            for flag in ("ready_set.actions.hydrogens=False", "ready_set.actions.ligands=True",
                         "ready_set.actions.optimise_ligand_geometry=False",
                         f"ready_set.input.cif_file_name={run / 'ReadySet' / 'input_dictionaries.cif'}"):
                self.assertIn(flag, readyset_args)
            self.assertEqual(result.readyset.curated_codes, ("DU",))
            self.assertEqual(result.readyset.dictionaries, (run / "ReadySet" / "DU.cif",))
            self.assertIn(str(run / "ReadySet" / "DU.cif"), received)
            self.assertEqual((run / "readyset_model.pdb").read_bytes(),
                             (run / "start_model.pdb").read_bytes())
            self.assertIn(str(mtz.resolve()), received)
            self.assertIn("xray_data.r_free_flags.generate=False", received)

            site = result.round.sites[0]
            self.assertTrue(site.calibrated)
            self.assertEqual((site.site, site.atom_name, site.element), ("B:101", "AG", "AG"))
            self.assertEqual((site.occupancy, site.b_factor), (0.9, 50.0))  # B from step 1
            self.assertEqual(site.f_prime, -0.12)
            self.assertEqual(site.f_double_prime_refined, 2.6)
            self.assertEqual(site.f_double_prime_theory, 2.03)
            self.assertAlmostEqual(site.deviation, 0.57)
            self.assertTrue(site.within_tolerance)
            self.assertEqual(result.round.warnings, ())
            self.assertEqual(result.round.map_coefficients.name, "refined_001_map_coeffs.mtz")

            report = json.loads(result.report_path.read_text())
            self.assertEqual((report["workflow"], report["schema_version"]), ("autofdp", 2))
            self.assertEqual(report["status"], STATUS_CONVERGED)
            self.assertEqual(report["steps"]["occupancy_adjustment"], "done")
            self.assertEqual(report["occupancy"]["sites"][0]["reason"],
                             "within +/-1 e- of theory without adjustment")
            self.assertEqual(report["final_model"]["relative_path"], "final_model.pdb")
            self.assertEqual((run / "final_model.pdb").read_bytes(),
                             (run / "readyset_model.pdb").read_bytes())
            # Inputs are anchored to the dataset folder (portable); the
            # absolute path is provenance only.
            self.assertEqual(
                {key: report["inputs"]["model"][key] for key in ("anchor", "relative_path", "sha256")},
                {"anchor": "dataset",
                 "relative_path": "AutoMR/run_001/AutoRefine/round_002/refined_001.pdb",
                 "sha256": file_sha256(model)},
            )
            self.assertEqual(report["inputs"]["model"]["absolute_path_at_creation"], str(model.resolve()))
            self.assertEqual(
                (report["inputs"]["mtz"]["relative_path"], report["inputs"]["mtz"]["sha256"]),
                ("staraniso_alldata-unique.mtz", file_sha256(mtz)),
            )
            self.assertEqual(report["reflections"]["wavelength"], 1.00743)
            self.assertEqual(
                (report["b_start"]["start_model"]["anchor"], report["b_start"]["start_model"]["relative_path"]),
                ("run", "start_model.pdb"),
            )
            self.assertEqual(report["b_start"]["sites"][0]["neighbour_mean_b"], 50.0)
            entry = report["rounds"][0]
            self.assertEqual(entry["directory"], "round_001")
            self.assertEqual(entry["strategy"], ["group_anomalous"])
            self.assertEqual(entry["input_model"]["relative_path"], "readyset_model.pdb")
            self.assertIsNone(entry["adjusted_site"])
            self.assertEqual(entry["status"], STATUS_ROUND_DONE)
            self.assertEqual(entry["starting_groups"][SILVER], {"f_prime": -0.12, "f_double_prime": 0.0})
            self.assertEqual(entry["sites"][0]["deviation"], 0.57)
            self.assertEqual(entry["held_parameter_changes"]["max_b_change"], 0.0)
            self.assertEqual(entry["refined_model"]["relative_path"], "round_001/refined_001.pdb")
            self.assertEqual(report["steps"]["readyset"], "done")
            self.assertEqual(report["readyset"]["mode"], "updated model")
            self.assertEqual(report["readyset"]["curated_residues"], ["DU"])
            self.assertEqual(report["readyset"]["model"]["relative_path"], "readyset_model.pdb")
            self.assertEqual(
                [item["relative_path"] for item in report["readyset"]["refinement_dictionaries"]],
                ["ReadySet/DU.cif"],
            )
            self.assertEqual(report["readyset"]["model_changes"]["max_b_change"], 0.0)
            log = result.log_path.read_text()
            for text in ("Step 1 (B start)", "B 379.99 -> 50.00", "ReadySet (no hydrogens)",
                         "Step 2", "deviation +0.57"):
                self.assertIn(text, log)

            # --max-iterations 0: round_001 only; a site outside tolerance is for review.
            again = _execute(_plan(root, model, mtz), root, refine=make_fdp_refine(root, {SILVER: 5.0}),
                             max_iterations=0)
            self.assertEqual(again.run_directory.name, "run_002")
            self.assertEqual(len(again.rounds), 1)
            self.assertAlmostEqual(again.round.sites[0].deviation, 2.97)
            self.assertFalse(again.round.sites[0].within_tolerance)
            self.assertIn("1 of 1 site(s) more than 1 e- from theory", again.round.message)
            self.assertEqual(again.status, STATUS_REVIEW)
            self.assertEqual(again.message, "B:101 AG: not adjusted (--max-iterations 0)")

    def test_record_resolves_after_the_dataset_folder_moves(self):
        from nasolve.run_context import resolve_artifact_path

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            result = _execute(_plan(root, model, mtz), root)
            moved = root / "moved" / "7CU_Ag"
            moved.parent.mkdir()
            mtz.parent.rename(moved)
            self.assertFalse(model.exists())
            run = moved / "AutoFDP" / "run_001"
            report = json.loads((run / "report.json").read_text())
            references = [report["inputs"]["model"], report["inputs"]["mtz"],
                          report["b_start"]["start_model"], report["rounds"][0]["refined_model"]]
            for reference in references:
                with self.subTest(reference=reference["relative_path"]):
                    path = resolve_artifact_path(reference, run)
                    self.assertIsNotNone(path)
                    self.assertTrue(str(path).startswith(str(moved.resolve())))
                    self.assertEqual(file_sha256(path), reference["sha256"])

    def test_theory_override_tolerance_and_output_folder(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            plan = _plan(root, model, mtz, theory_overrides={"AG": 4.5}, output_root=root / "elsewhere")
            result = _execute(plan, root, refine=make_fdp_refine(root, {SILVER: 5.0}),
                              tolerance=0.4, cycles=2)
            self.assertEqual(result.run_directory, (root / "elsewhere" / "run_001").resolve())
            site = result.rounds[0].sites[0]
            self.assertEqual((site.f_double_prime_theory, site.theory_source), (4.5, "user (--theory)"))
            self.assertAlmostEqual(site.deviation, 0.5)
            self.assertFalse(site.within_tolerance)
            self.assertIn("number_of_macro_cycles = 2", result.round.parameters.read_text())
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["theory"]["AG"]["f_double_prime_henke"], 2.03)
            # f' stays the calculated value even when f'' theory is overridden.
            self.assertEqual(report["rounds"][0]["starting_groups"][SILVER]["f_prime"], -0.12)
            with self.assertRaisesRegex(AutoFDPError, "HG, which is not an anomalous scatterer"):
                _plan(root, model, mtz, theory_overrides={"HG": 10.17})

    def test_every_metal_is_refined_but_site_limits_the_b_start(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)
            plan = _plan(root, model, mtz, sites=["B:101"])
            self.assertEqual(plan.selections, (SILVER, MERCURY))
            result = _execute(plan, root, sites=["B:101"])
            silver, mercury = result.round.sites
            self.assertTrue(silver.calibrated)
            self.assertFalse(mercury.calibrated)
            self.assertEqual((mercury.b_factor, mercury.f_prime, mercury.f_double_prime_theory),
                             (150.0, -12.5, 10.1))   # B unchanged, own element theory
            self.assertAlmostEqual(mercury.deviation, -1.1)
            self.assertIn(f'selection = "{MERCURY}"', result.round.parameters.read_text())

    def test_missing_refined_fdp_or_changed_held_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            result = _execute(_plan(root, model, mtz), root, refine=make_fdp_refine(root, {}))
            self.assertEqual(result.status, STATUS_FAILED)
            self.assertIn("no refined f''", result.message)
            self.assertIsNone(result.round.sites[0].deviation)
            report = json.loads(result.report_path.read_text())
            self.assertEqual((report["status"], report["rounds"][0]["status"]), (STATUS_FAILED, STATUS_FAILED))
            self.assertEqual(report["steps"]["occupancy_adjustment"], "not run")
            self.assertNotIn("final_model", report)
            self.assertIsNone(result.final_model)

            changed = _execute(_plan(root, model, mtz), root, refine=make_fdp_refine(root, change_b=True))
            self.assertIn("Held parameters changed", changed.round.warnings[0])
            self.assertEqual(
                held_parameter_changes(changed.readyset.model, changed.round.model_path)["max_b_change"],
                5.0,
            )


class ReadySetStepTests(unittest.TestCase):
    """ReadySet (no hydrogens) makes the restraint dictionaries for the f''
    round from the model alone; it must not change the model."""

    def test_generated_and_user_dictionaries(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            # A --cif for a residue ReadySet regenerates, and one it does not.
            regenerated = root / "pdu.cif"
            regenerated.write_text(
                "data_comp_PDU\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n"
                "_chem_comp_atom.type_symbol\nPDU C1 C\n"
            )
            extra = root / "xyz.cif"
            extra.write_text(
                "data_comp_XYZ\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n"
                "_chem_comp_atom.type_symbol\nXYZ C1 C\n"
            )
            plan = _plan(root, model, mtz, restraints=[regenerated, extra])
            result = _execute(plan, root, ready_set=make_ready_set(root, generated=("PDU",)))
            run = result.run_directory
            names = [path.name for path in result.readyset.dictionaries]
            # ReadySet's start_model.ligands.cif goes on as readyset_supplement.cif
            # (PostMR's precedence: curated overrides filtered out).
            self.assertEqual(names, ["readyset_supplement.cif", "DU.cif", "user_xyz.cif"])
            self.assertEqual(result.readyset.generated,
                             run / "ReadySet" / "start_model.ligands.cif")
            received = json.loads((run / "round_001" / "received_args.json").read_text())
            self.assertIn(str(run / "ReadySet" / "readyset_supplement.cif"), received)
            self.assertNotIn(str(run / "ReadySet" / "user_pdu.cif"), received)  # not twice
            from nasolve.autofdp import _component_codes

            self.assertEqual(_component_codes(run / "ReadySet" / "input_dictionaries.cif"),
                             {"DU", "PDU", "XYZ"})

    def test_no_unknown_residues_keeps_the_start_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            result = _execute(_plan(root, model, mtz), root, ready_set=make_ready_set(root, mode="noop"))
            self.assertEqual(result.status, STATUS_CONVERGED)
            self.assertEqual(result.readyset.mode, "no unknown residues")
            self.assertIsNone(result.readyset.generated)
            self.assertEqual((result.run_directory / "readyset_model.pdb").read_bytes(),
                             result.b_start.start_model.read_bytes())

    def test_readyset_problems_stop_before_refinement(self):
        cases = (("fail", "ReadySet failed with status 1"),
                 ("change_b", "ReadySet changed the model"),
                 ("add_h", "ReadySet added 1 hydrogen"))
        for mode, message in cases:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                model, mtz = make_inputs(root)
                with self.assertRaisesRegex(AutoFDPError, message):
                    _execute(_plan(root, model, mtz), root, ready_set=make_ready_set(root, mode=mode))
                run = mtz.parent.resolve() / "AutoFDP" / "run_001"
                report = json.loads((run / "report.json").read_text())
                self.assertEqual((report["status"], report["steps"]["readyset"]),
                                 (STATUS_FAILED, "failed"))
                self.assertFalse((run / "round_001").exists())
                self.assertIn("Failed:", (run / "autofdp.log").read_text())


class AutoFDPCommandTests(unittest.TestCase):
    def _main(self, root: Path, arguments: list[str], dump=None, refine=None):
        """Run ``nasolve autofdp``; the last value is a mock that fails the
        test if AutoFDP looks for Coot (it never opens Coot)."""
        installation = PhenixInstallation(
            root=None, setup_script=None, version=PHENIX_21,
            executables={"phenix.refine": refine or make_fdp_refine(root),
                         "phenix.ready_set": make_ready_set(root),
                         "phenix.mtz.dump": dump or make_mtz_dump(root, wavelength=1.00743)},
            environment=ENV, source="test",
        )
        out, err = io.StringIO(), io.StringIO()
        with (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli.discover_phenix", return_value=installation),
            patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
            patch("nasolve.cli.discover_coot",
                  side_effect=AssertionError("AutoFDP must not look for Coot")) as coot_lookup,
            contextlib.redirect_stdout(out), contextlib.redirect_stderr(err),
        ):
            code = main(["autofdp", *arguments])
        return code, out.getvalue(), err.getvalue(), coot_lookup

    def test_command_runs_both_steps_without_coot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            code, out, err, coot_lookup = self._main(root, [str(model), str(mtz)])
            self.assertEqual(code, 0, err)
            for text in (
                "Observations: F(+),SIGF(+),F(-),SIGF(-); Free-R FreeR_flag (test value 0)",
                "Wavelength: 1.00743 A (MTZ header (phenix.mtz.dump))",
                "AG f' -0.12, f'' 2.03 (cctbx-henke)",
                "B:101 AG (AG), occupancy 0.90: B 379.99 -> 50.00",
                "f''-only round (round_001, 3 cycles",
                "Theory: AG cctbx-henke",
                "AUTOFDP_CONVERGED",
                "every site already within +/-1 e- of theory; no occupancy changed",
                "Model with the final occupancies: ",
                "ReadySet (no hydrogens): updated model; atoms, coordinates and B unchanged",
                "curated dictionaries given: DU",
                "restraints for phenix.refine: ReadySet/DU.cif",
                "Not implemented yet: the identity check, B confirmation",
            ):
                self.assertIn(text, out)
            row = next(line for line in out.splitlines() if line.startswith("B:101 AG"))
            for value in ("0.90", "50.00", "-0.12", "2.60", "2.03", "+0.57", "within +/-1"):
                self.assertIn(value, row)
            # The maps of an f''-only round are not diagnostic: no Coot at all.
            run = mtz.parent.resolve() / "AutoFDP" / "run_001"
            coot_lookup.assert_not_called()
            self.assertNotIn("Coot", out)
            self.assertFalse((run / "CootGUI").exists())

            # Either order, with a theory override: a new run folder. The old
            # --no-coot is still accepted and changes nothing.
            code, out, err, coot_lookup = self._main(
                root, [str(mtz), str(model), "--no-coot", "--theory", "AG=2.06"]
            )
            self.assertEqual(code, 0, err)
            coot_lookup.assert_not_called()
            self.assertIn("Theory: AG user (--theory)", out)
            self.assertIn(
                "Model with the final occupancies: "
                + str(mtz.parent.resolve() / "AutoFDP" / "run_002" / "final_model.pdb"), out)
            with contextlib.redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit) as exit:
                main(["autofdp", str(model), str(mtz), "--coot", "/opt/coot"])
            self.assertEqual(exit.exception.code, 2)
            self.assertIn("unrecognized arguments: --coot", err.getvalue())

    def test_plan_refines_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            code, out, err, coot_lookup = self._main(root, [str(model), str(mtz), "--plan"])
            self.assertEqual(code, 0, err)
            self.assertIn("B 379.99 -> 50.00", out)
            self.assertIn("Then: ReadySet on the start model (no hydrogens)", out)
            self.assertIn("one f''-only phenix.refine round of 3 cycles", out)
            self.assertIn("then up to 5 occupancy adjustment(s) per site", out)
            self.assertFalse((root / "readyset_args.json").exists())
            self.assertIn("Plan only", out)
            coot_lookup.assert_not_called()
            self.assertFalse((mtz.parent / "AutoFDP").exists())
            self.assertFalse((root / "received_args.json").exists())

    def test_command_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            cases = (
                ([str(model), str(model)], 1, "Give one model (.pdb) and one reflection file (.mtz)"),
                ([str(model), str(mtz), "--theory", "AG"], 1, "--theory must look like HG=10.17"),
                ([str(model), str(mtz), "--cycles", "0"], 1, "--cycles must be at least 1"),
                ([str(model), str(mtz), "--cif", str(root / "none.cif")], 1, "No restraint file"),
                ([str(model), str(mtz), "--element", "B:101=HG"], 2, "--element is not implemented yet"),
                ([str(model), str(mtz), "--max-iterations", "-1"], 1, "--max-iterations must be 0 or more"),
                ([str(model), str(mtz), "--occupancy-floor", "1"], 1,
                 "--occupancy-floor must be above 0 and below 1"),
            )
            for arguments, expected, message in cases:
                with self.subTest(arguments=arguments[2:]):
                    code, _, err, _ = self._main(root, arguments)
                    self.assertEqual(code, expected)
                    self.assertIn(message, err)
            self.assertFalse((mtz.parent / "AutoFDP").exists())
            code, _, err, _ = self._main(root, [str(model), str(mtz)], dump=make_mtz_dump(root))
            self.assertEqual(code, 1)
            self.assertIn("give the wavelength with --wavelength", err)


if __name__ == "__main__":
    unittest.main()
