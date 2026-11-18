"""AutoFDP's final round (--final-eff): after convergence, a phenix.refine
round with every setting of a manual .eff, AutoFDP's model, the metal sites'
occupancy held and f'/f'' fixed at the calibrated values, opened in Coot.

phenix.refine, phenix.ready_set, phenix.mtz.dump, phenix.python and Coot are
simulated. The stand-in phenix.refine acts as an f'' round until it is given
final_input.eff, then as the final refinement.
"""

import contextlib
import difflib
import io
import json
import stat
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.autofdp import (
    STATUS_CONVERGED,
    STATUS_FAILED,
    STATUS_FINAL_DONE,
    STATUS_REVIEW,
    AutoFDPError,
    FdpSiteResult,
    read_final_eff,
    write_final_inputs,
)
from nasolve.cli import main
from nasolve.config import AppConfig
from nasolve.coot_runtime import CootInstallation
from nasolve.coot_view import CootViewError, launch_coot_files
from nasolve.phenix_runtime import PhenixInstallation
from nasolve.phil_text import scan, starred, unquote

from .test_autofdp_round import (
    ENV,
    SILVER,
    _execute,
    _plan,
    make_inputs,
    make_ready_set,
)
from .test_autorefine import PHENIX_21, make_mtz_dump


def write_manual_eff(root: Path, model: Path, mtz: Path, *, reflection: Path | None = None,
                     extra_model: str = "", self_reference: bool = False) -> Path:
    """A GUI-style .eff of a manual refinement (like refine_77.eff)."""
    restraints = root / "restraints"
    restraints.mkdir(exist_ok=True)
    (restraints / "DU.cif").write_text("data_comp_DU\n")
    (restraints / "pairs.phil").write_text("refinement { }\n")
    eff = root / "project_data" / "refine_12.eff"
    eff.parent.mkdir(exist_ok=True)
    phil_self = f'  phil_files = "{eff}"\n' if self_reference else ""
    eff.write_text(
        "data_manager {\n"
        f'  restraint_files = "{restraints / "DU.cif"}"\n'
        f'  phil_files = "{restraints / "pairs.phil"}"\n'
        + phil_self +
        "  default_phil = None\n"
        "  model {\n"
        f'    file = "{model}"\n'
        "    type = *x_ray neutron electron reference\n"
        "  }\n"
        + extra_model +
        f'  default_model = "{model}"\n'
        "  miller_array {\n"
        f'    file = "{reflection or mtz}"\n'
        "    labels {\n"
        '      name = "F(+),SIGF(+),F(-),SIGF(-)"\n'
        "    }\n"
        '    user_selected_labels = "F(+),SIGF(+),F(-),SIGF(-)"\n'
        "  }\n"
        "}\n"
        "refinement {\n"
        "  refine {\n"
        "    strategy = *individual_sites *individual_sites_real_space rigid_body \\\n"
        "               *individual_adp group_adp tls *occupancies *group_anomalous \\\n"
        "               den\n"
        "    occupancies {\n"
        "      individual = None\n"
        "      remove_selection = chain A resid 12\n"
        "    }\n"
        "    anomalous_scatterers {\n"
        "      group {\n"
        "        selection = chain B resid 101\n"
        "        f_prime = 0\n"
        "        f_double_prime = 0\n"
        "        refine = *f_prime *f_double_prime\n"
        "      }\n"
        "    }\n"
        "  }\n"
        "  main {\n"
        "    number_of_macro_cycles = 5\n"
        "    random_seed = 2679941\n"
        "  }\n"
        "  gui {\n"
        "    base_output_dir = /data/Structural \\\n"
        "                      Biology/project\n"
        "  }\n"
        "}\n"
        "output {\n"
        '  prefix = """manual_refine"""\n'
        "  serial = 12\n"
        "  overwrite = False\n"
        "}\n"
    )
    return eff


def make_final_refine(root: Path, fdp: float = 2.6, *, final: str = "ok") -> Path:
    """Stand-in phenix.refine. f'' rounds: echo each anomalous group with f''
    ``fdp`` (Ag only). Final round (given final_input.eff): ``ok`` writes
    final_001.pdb/.mtz (B of the Ag up by 5) and echoes the fixed groups;
    ``fail`` exits 1; ``refuse`` fails the dry run; ``drift`` also moves the
    Ag occupancy and reports another f''."""
    executable = root / "phenix.refine"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, re, shutil, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "is_final = any(a.endswith('final_input.eff') for a in args)\n"
        f"mode = {final!r}\n"
        "if is_final:\n"
        "    eff = next(Path(a) for a in args if a.endswith('final_input.eff')).read_text()\n"
        "    loaded = [Path(m) for m in re.findall(r'phil_files = \"(.*?)\"', eff)]\n"
        "    texts = [eff, *(path.read_text() for path in loaded)]\n"
        "    given = next(Path(a) for a in args if a.endswith('.pdb')).resolve()\n"
        "    named = {Path(v).resolve() for t in texts for v in re.findall(r'= \"([^\"]+[.]pdb)\"', t)}\n"
        "    if sum(text.count('  model {') for text in texts) > 0 or any(\n"
        "            path.is_file() and path != given for path in named):\n"
        "        print('Sorry: Wrong number of models of each type supplied.')\n"
        "        raise SystemExit(1)\n"
        "    if any('anomalous_scatterers' in text for text in texts[1:]):\n"
        "        print('Sorry: an atom is in more than one anomalous group')\n"
        "        raise SystemExit(1)\n"
        "if '--dry_run' in args:\n"
        "    if is_final and mode == 'refuse':\n"
        "        print('Sorry: invalid selection')\n"
        "        raise SystemExit(1)\n"
        "    print('Dry-run parameters accepted')\n"
        "    raise SystemExit(0)\n"
        "Path('received_args.json').write_text(json.dumps(args))\n"
        "params = next(Path(a) for a in args if a.endswith('.params')).read_text()\n"
        "model = next(Path(a) for a in args if a.endswith('.pdb'))\n"
        "def echo(fdp_of):\n"
        "    for block in re.findall(r'group \\{(.*?)\\}', params, re.S):\n"
        "        selection = re.search(r'selection = \"(.*?)\"', block).group(1)\n"
        "        print('Anomalous scatterer group:')\n"
        "        print(f'  Selection: \"{selection}\"')\n"
        "        print('  Number of selected scatterers: 1')\n"
        "        print('  f_prime: ' + re.search(r'f_prime = (\\S+)', block).group(1))\n"
        "        print(f'  f_double_prime: {fdp_of(block)}')\n"
        "if not is_final:\n"
        "    shutil.copyfile(model, 'refined_001.pdb')\n"
        "    Path('refined_001_map_coeffs.mtz').write_bytes(b'maps')\n"
        "    print('Final: r_work = 0.2161 r_free = 0.2121 bonds = 0.012 angles = 1.82')\n"
        f"    echo(lambda block: {fdp!r})\n"
        "    raise SystemExit(0)\n"
        "if mode == 'fail':\n"
        "    print('Sorry: fixture failure')\n"
        "    raise SystemExit(1)\n"
        "lines = model.read_text().splitlines(keepends=True)\n"
        "for index, line in enumerate(lines):\n"
        "    if line.startswith('HETATM') and line[12:16].strip() == 'AG':\n"
        "        occupancy = float(line[54:60]) - (0.2 if mode == 'drift' else 0)\n"
        "        line = line[:54] + '%6.2f' % occupancy + '%6.2f' % (float(line[60:66]) + 5) + line[66:]\n"
        "        lines[index] = line\n"
        "Path('final_001.pdb').write_text(''.join(lines))\n"
        "Path('final_001.mtz').write_bytes(b'data and maps')\n"
        "print('Final: r_work = 0.1903 r_free = 0.2310 bonds = 0.008 angles = 1.10')\n"
        "echo(lambda block: 9.99 if mode == 'drift' else re.search(r'f_double_prime = (\\S+)', block).group(1))\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    python = root / "phenix.python"
    python.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' 'NASOLVE_SCATTERING_JSON:{\"AG\": [-0.12, 2.03], \"HG\": [-12.5, 10.1]}'\n"
    )
    python.chmod(python.stat().st_mode | stat.S_IXUSR)
    return executable


class ReadFinalEffTests(unittest.TestCase):
    def test_reads_what_the_final_round_keeps_and_replaces(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            manual = read_final_eff(write_manual_eff(root, model, mtz), mtz)
            self.assertEqual(manual.models, (str(model),))
            self.assertEqual(manual.reflection_files, (str(mtz),))
            self.assertEqual(manual.restraint_files, (str(root / "restraints" / "DU.cif"),))
            self.assertEqual(manual.phil_files, (str(root / "restraints" / "pairs.phil"),))
            self.assertEqual(manual.strategy, ("individual_sites", "individual_sites_real_space",
                                               "individual_adp", "occupancies", "group_anomalous"))
            self.assertEqual(manual.anomalous_groups, ("chain B resid 101",))
            self.assertEqual(manual.occupancy_exclusions, ("chain A resid 12",))
            self.assertEqual(manual.warnings, ())
            other = root / "other.mtz"
            other.write_bytes(b"x")
            warned = read_final_eff(write_manual_eff(root, model, mtz, reflection=other), mtz)
            self.assertIn("refines against", warned.warnings[0])

    def test_refuses_what_phenix_could_not_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            eff = write_manual_eff(root, model, mtz)
            plain = root / "plain.eff"
            for text, message in (
                ("data_manager {\n", "not a parameter file NASolve can read"),
                ('data_manager {\n  default_model = "x.pdb"\n}\n', "no refinement parameters"),
                ("refinement {\n  main {\n    number_of_macro_cycles = 3\n  }\n}\n",
                 "names no reflection file"),
            ):
                with self.subTest(text=text):
                    plain.write_text(text)
                    with self.assertRaisesRegex(AutoFDPError, message):
                        read_final_eff(plain, mtz)
            with self.assertRaisesRegex(AutoFDPError, "no such file"):
                read_final_eff(root / "none.eff", mtz)
            (root / "restraints" / "DU.cif").unlink()
            with self.assertRaisesRegex(AutoFDPError, "refers to files that do not exist: .*DU.cif"):
                read_final_eff(eff, mtz)


def silver_site(fdp: float | None = 2.6) -> FdpSiteResult:
    return FdpSiteResult(
        selection=SILVER, site="B:101", atom_name="AG", element="AG", calibrated=True,
        occupancy=0.9, b_factor=50.0, f_prime=-0.12, f_double_prime_start=0.0,
        f_double_prime_refined=fdp, f_double_prime_theory=2.03, theory_source="cctbx-henke",
        deviation=None if fdp is None else round(fdp - 2.03, 3), within_tolerance=True,
    )


class WriteFinalInputsTests(unittest.TestCase):
    def test_eff_copy_changes_only_model_groups_and_strategy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            reference = (
                "  model {\n"
                f'    file = "{root / "reference.pdb"}"\n'
                "    type = x_ray neutron electron *reference\n"
                "  }\n"
            )
            eff = write_manual_eff(root, model, mtz, extra_model=reference, self_reference=True)
            manual = read_final_eff(eff, mtz)
            written, parameters, edits = write_final_inputs(manual, root / "final_round", [silver_site()])
            text = written.read_text()
            self.assertTrue(text.startswith("# NASolve AutoFDP final round: copy of "))
            items = scan(text)
            paths = [item.path for item in items]
            self.assertNotIn("data_manager.default_model", paths)
            self.assertNotIn("refinement.refine.anomalous_scatterers", paths)
            models = [unquote(item.value) for item in items if item.path == "data_manager.model.file"]
            self.assertEqual(models, [str(root / "reference.pdb")])   # the reference model stays
            self.assertEqual([unquote(item.value) for item in items if item.path == "data_manager.phil_files"],
                             [str(root / "restraints" / "pairs.phil")])
            strategy = next(item for item in items if item.path == "refinement.refine.strategy")
            self.assertEqual(starred(strategy.value), ("individual_sites", "individual_sites_real_space",
                                                       "individual_adp", "occupancies"))
            self.assertIn("               *individual_adp group_adp tls *occupancies group_anomalous \\\n", text)
            self.assertIn("    number_of_macro_cycles = 5\n", text)
            self.assertIn("      remove_selection = chain A resid 12\n", text)
            self.assertIn('  prefix = """manual_refine"""\n', text)
            self.assertEqual(edits["removed_models"], [str(model)])
            self.assertEqual(edits["removed_anomalous_groups"], ["chain B resid 101"])
            self.assertEqual(edits["removed_phil_files"], [str(eff)])
            self.assertEqual(edits["strategy_after"], ["individual_sites", "individual_sites_real_space",
                                                       "individual_adp", "occupancies"])

            # The metal sites: occupancy held, f' and f'' fixed, outputs named final_001.
            extra = scan(parameters.read_text())
            values = {item.path: item.value for item in extra if item.kind == "assign"}
            self.assertEqual(unquote(values["refinement.refine.occupancies.remove_selection"]), SILVER)
            self.assertEqual(unquote(values["refinement.refine.anomalous_scatterers.group.selection"]), SILVER)
            self.assertEqual(values["refinement.refine.anomalous_scatterers.group.f_prime"], "-0.12")
            self.assertEqual(values["refinement.refine.anomalous_scatterers.group.f_double_prime"], "2.6")
            self.assertEqual(values["refinement.refine.anomalous_scatterers.group.refine"], "None")
            self.assertEqual((unquote(values["output.prefix"]), values["output.serial"],
                              values["output.overwrite"]), ("final", "1", "True"))

    def test_only_the_listed_lines_differ(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            eff = write_manual_eff(root, model, mtz)
            written, _, _ = write_final_inputs(read_final_eff(eff, mtz), root / "final_round", [silver_site()])
            changed = [line for line in difflib.ndiff(eff.read_text().splitlines(),
                                                       written.read_text().splitlines())
                       if line[:1] in "+-"]
            removed = [line[2:] for line in changed if line.startswith("- ")]
            added = [line[2:] for line in changed if line.startswith("+ ")]
            self.assertEqual(sorted(removed), sorted([
                "  model {", f'    file = "{model}"', "    type = *x_ray neutron electron reference",
                "  }", f'  default_model = "{model}"',
                "               *individual_adp group_adp tls *occupancies *group_anomalous \\",
                "    anomalous_scatterers {", "      group {", "        selection = chain B resid 101",
                "        f_prime = 0", "        f_double_prime = 0",
                "        refine = *f_prime *f_double_prime", "      }", "    }",
            ]))
            self.assertEqual([line for line in added if not line.startswith("#")],
                             ["               *individual_adp group_adp tls *occupancies group_anomalous \\"])


class FinalRoundTests(unittest.TestCase):
    def test_converged_run_ends_with_the_final_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_final_refine(root)
            eff = write_manual_eff(root, model, mtz)
            result = _execute(_plan(root, model, mtz, refine=refine, final_eff=eff), root, refine=refine)
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)
            final = result.final
            self.assertEqual(final.status, STATUS_FINAL_DONE, final.message)
            run = result.run_directory
            self.assertEqual(final.directory, run / "final_round")
            self.assertEqual(list(final.command)[1:], [
                str(run / "final_model.pdb"), str(run / "final_round" / "final_input.eff"),
                str(run / "final_round" / "autofdp_final.params"),
            ])
            self.assertEqual(json.loads((run / "final_round" / "received_args.json").read_text()),
                             list(final.command)[1:])
            self.assertEqual((final.r_work, final.r_free), (0.1903, 0.231))
            self.assertEqual(final.model_path, run / "final_round" / "final_001.pdb")
            self.assertEqual(final.map_coefficients, run / "final_round" / "final_001.mtz")
            site = final.sites[0]
            self.assertEqual(
                (site.occupancy, site.occupancy_after, site.f_prime, site.f_double_prime,
                 site.f_double_prime_reported, site.b_before, site.b_after),
                (0.9, 0.9, -0.12, 2.6, 2.6, 50.0, 55.0),
            )
            self.assertEqual(final.warnings, ())
            self.assertIn("final round done with refine_12.eff's settings", result.message)

            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["status"], STATUS_CONVERGED)
            self.assertEqual(report["steps"]["final_round"], "done")
            self.assertEqual(report["settings"]["final_eff"], str(eff))
            entry = report["final_round"]
            self.assertEqual(entry["status"], STATUS_FINAL_DONE)
            self.assertEqual(entry["eff_written"]["relative_path"], "final_round/final_input.eff")
            self.assertEqual(entry["refined_model"]["relative_path"], "final_round/final_001.pdb")
            self.assertEqual(entry["edits"]["removed_anomalous_groups"], ["chain B resid 101"])
            self.assertEqual(entry["sites"][0]["f_double_prime"], 2.6)
            self.assertIn("sha256", entry["eff_source"])
            self.assertIn("Final round (final_round, settings of", result.log_path.read_text())

    def test_job_output_eff_that_names_the_gui_eff(self):
        # Refine_85/PdUPdU-Hg_refine_085.eff (written by phenix.refine) names the
        # GUI's .phenix/project_data/refine_85.eff in phil_files and default_phil;
        # that file has the model, the anomalous groups and the output name.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            gui_eff = write_manual_eff(root, model, mtz)
            gui_copy = root / "project_data" / "refine_85.eff"
            gui_copy.write_text(gui_eff.read_text())
            anomalous_phil = root / "restraints" / "narestraints_anomalous.phil"
            anomalous_phil.write_text(
                "refinement {\n  refine {\n"
                "    strategy = *individual_sites *individual_adp *group_anomalous\n"
                "    anomalous_scatterers {\n      group {\n        selection = chain B resid 101\n"
                "      }\n    }\n  }\n"
                "  pdb_interpretation {\n    secondary_structure {\n      enabled = True\n    }\n  }\n}\n"
            )
            text = gui_eff.read_text().replace(
                "  default_phil = None\n",
                f'  phil_files = "{gui_copy}"\n  phil_files = "{anomalous_phil}"\n'
                f'  default_phil = "{gui_copy}"\n',
            ).replace(
                "                      Biology/project\n",
                "                      Biology/project\n"
                "    migration {\n      refinement {\n        input {\n          pdb {\n"
                f'            file_name = "{model}"\n'
                "          }\n          xray_data {\n"
                f'            file_name = "{mtz}"\n'
                "          }\n        }\n      }\n    }\n",
            )
            output_eff = root / "Refine_85" / "PdUPdU-Hg_refine_085.eff"
            output_eff.parent.mkdir()
            output_eff.write_text(text)
            manual = read_final_eff(output_eff, mtz)
            self.assertEqual(manual.job_phil_files, (str(gui_copy),))
            plan_out = io.StringIO()
            with contextlib.redirect_stdout(plan_out):
                from nasolve.cli import _print_autofdp_plan
                _print_autofdp_plan(_plan(root, model, mtz, refine=make_final_refine(root),
                                          final_eff=output_eff))
            self.assertIn(f"left out: {gui_copy} (a refinement job's own parameters, already in the .eff)",
                          plan_out.getvalue())
            refine = make_final_refine(root)
            result = _execute(_plan(root, model, mtz, refine=refine, final_eff=output_eff), root,
                              refine=refine)
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)
            self.assertEqual(result.final.status, STATUS_FINAL_DONE)
            written = (result.run_directory / "final_round" / "final_input.eff").read_text()
            names = [unquote(item.value) for item in scan(written) if item.path == "data_manager.phil_files"]
            copy = result.run_directory / "final_round" / "phil_files" / "03_narestraints_anomalous.phil"
            self.assertEqual(names, [str(root / "restraints" / "pairs.phil"), str(copy)])
            self.assertNotIn("default_phil", [item.path.split(".")[-1] for item in scan(written)])
            copied = copy.read_text()
            self.assertNotIn("anomalous_scatterers", copied)
            self.assertIn("*individual_adp group_anomalous\n", copied)
            self.assertIn("enabled = True", copied)
            edits = result.final.edits
            self.assertEqual(edits["removed_phil_files"], [str(gui_copy)])
            self.assertTrue(edits["removed_default_phil"])
            # The GUI's legacy record of the old model is gone; its MTZ entry stays.
            self.assertEqual(edits["removed_model_references"],
                             ["refinement.gui.migration.refinement.input.pdb.file_name"])
            self.assertNotIn(str(model), written)
            self.assertIn(f'            file_name = "{mtz}"\n', written)
            self.assertEqual(edits["copied_phil_files"], [{
                "original": str(anomalous_phil), "copy": str(copy),
                "removed": ["group_anomalous in strategy", "refinement.refine.anomalous_scatterers"],
            }])
            self.assertEqual(anomalous_phil.read_text().count("anomalous_scatterers"), 1)  # untouched

    def test_no_final_round_without_convergence_or_without_an_eff(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_final_refine(root, fdp=5.0)
            eff = write_manual_eff(root, model, mtz)
            result = _execute(_plan(root, model, mtz, refine=refine, final_eff=eff), root,
                              refine=refine, max_iterations=0)
            self.assertEqual(result.status, STATUS_REVIEW)
            self.assertIsNone(result.final)
            self.assertFalse((result.run_directory / "final_round").exists())
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["steps"]["final_round"], "skipped: not converged")
            plain = _execute(_plan(root, model, mtz, refine=make_final_refine(root)), root,
                             refine=make_final_refine(root))
            self.assertIsNone(plain.final)
            self.assertEqual(json.loads(plain.report_path.read_text())["steps"]["final_round"],
                             "not requested")

    def test_final_round_problems(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            eff = write_manual_eff(root, model, mtz)
            for mode, log_name in (("fail", "phenix.refine.log"), ("refuse", "phenix.refine.preflight.log")):
                with self.subTest(mode=mode):
                    refine = make_final_refine(root, final=mode)
                    result = _execute(_plan(root, model, mtz, refine=refine, final_eff=eff), root,
                                      refine=refine)
                    self.assertEqual(result.status, STATUS_FAILED)
                    self.assertTrue(result.message.startswith(
                        "every site converged, but phenix.refine failed in the final round"))
                    self.assertTrue((result.run_directory / "final_round" / log_name).is_file())
                    report = json.loads(result.report_path.read_text())
                    self.assertEqual((report["status"], report["steps"]["final_round"]),
                                     (STATUS_FAILED, "failed"))
                    self.assertTrue((result.run_directory / "final_model.pdb").is_file())
            refine = make_final_refine(root, final="drift")
            drifted = _execute(_plan(root, model, mtz, refine=refine, final_eff=eff), root, refine=refine)
            self.assertEqual(drifted.final.status, STATUS_FINAL_DONE)
            self.assertEqual(drifted.final.warnings, (
                "B:101 AG occupancy moved from 0.90 to 0.70",
                "phenix.refine's log shows f'' 9.99 for B:101 AG, not the 2.60 set",
            ))


class FinalCommandTests(unittest.TestCase):
    def _main(self, root: Path, arguments: list[str], refine: Path, launch=None):
        installation = PhenixInstallation(
            root=None, setup_script=None, version=PHENIX_21,
            executables={"phenix.refine": refine, "phenix.ready_set": make_ready_set(root),
                         "phenix.mtz.dump": make_mtz_dump(root, wavelength=1.00743)},
            environment=ENV, source="test",
        )
        coot = root / "coot"
        coot.write_text("#!/bin/sh\n")
        launch = launch or (lambda model, maps, *args, **kwargs: SimpleNamespace(
            pid=4321, model_path=model, map_path=maps, log_path=Path("coot_gui.log"),
        ))
        out, err = io.StringIO(), io.StringIO()
        with (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli.discover_phenix", return_value=installation),
            patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
            patch("nasolve.cli.discover_coot", return_value=CootInstallation(coot, "test", "fixture")),
            patch("nasolve.cli.remember_coot"),
            patch("nasolve.cli.launch_coot_files", side_effect=launch) as launched,
            contextlib.redirect_stdout(out), contextlib.redirect_stderr(err),
        ):
            code = main(["autofdp", *arguments])
        return code, out.getvalue(), err.getvalue(), launched

    def test_final_round_is_printed_and_opened_in_coot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            eff = write_manual_eff(root, model, mtz)
            code, out, err, launched = self._main(
                root, [str(model), str(mtz), "--final-eff", str(eff)], make_final_refine(root)
            )
            self.assertEqual(code, 0, err)
            run = mtz.parent.resolve() / "AutoFDP" / "run_001"
            for text in (
                f"Final round settings (--final-eff): {eff}",
                "(group_anomalous will be unselected)",
                "its anomalous groups (chain B resid 101) -> replaced by the metal sites",
                "occupancies held as in the .eff: chain A resid 12",
                "AutoFDP final round (final_round, settings of refine_12.eff; metal occupancy held, "
                "f' and f'' fixed)",
                "R-work 0.1903, R-free 0.231",
                f"Final model: {run / 'final_round' / 'final_001.pdb'}",
                "Opened Coot (PID 4321) with final_round/final_001.pdb and final_001.mtz",
            ):
                self.assertIn(text, out)
            self.assertLess(out.index("AutoFDP final round ("), out.index("AUTOFDP_CONVERGED"))
            row = next(line for line in out.splitlines() if line.startswith("B:101 AG") and "55.00" in line)
            for value in ("0.90", "-0.12", "2.60", "50.00", "55.00"):
                self.assertIn(value, row)
            args, kwargs = launched.call_args
            self.assertEqual(args[:2], (run / "final_round" / "final_001.pdb", run / "final_round" / "final_001.mtz"))
            self.assertEqual(kwargs["working_directory"], run / "CootGUI")
            self.assertEqual(kwargs["dictionaries"], [root / "restraints" / "DU.cif"])

    def test_no_coot_coot_problems_plan_and_not_converged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            eff = write_manual_eff(root, model, mtz)
            code, out, err, launched = self._main(
                root, [str(model), str(mtz), "--final-eff", str(eff), "--no-coot"], make_final_refine(root)
            )
            self.assertEqual(code, 0, err)
            launched.assert_not_called()
            self.assertIn("Maps: ", out)

            def broken(*args, **kwargs):
                raise CootViewError("fixture: no display")
            code, out, err, _ = self._main(
                root, [str(model), str(mtz), "--final-eff", str(eff)], make_final_refine(root), launch=broken
            )
            self.assertEqual(code, 0)
            self.assertIn("Coot view error: fixture: no display", err)

            code, out, err, launched = self._main(
                root, [str(model), str(mtz), "--final-eff", str(eff), "--plan"], make_final_refine(root)
            )
            self.assertEqual(code, 0, err)
            self.assertIn("If every site converges: a final round with the settings of refine_12.eff", out)
            launched.assert_not_called()

            code, out, err, launched = self._main(
                root, [str(model), str(mtz), "--final-eff", str(eff), "--max-iterations", "0"],
                make_final_refine(root, fdp=5.0),
            )
            self.assertEqual(code, 2, err)
            self.assertIn("AutoFDP final round: not run, because not every site converged", out)
            launched.assert_not_called()

            code, _, err, _ = self._main(
                root, [str(model), str(mtz), "--coot", "/opt/coot", "--no-coot"], make_final_refine(root)
            )
            self.assertEqual(code, 1)
            self.assertIn("--coot and --no-coot cannot be combined", err)
            code, _, err, _ = self._main(
                root, [str(model), str(mtz), "--final-eff", str(root / "none.eff")], make_final_refine(root)
            )
            self.assertEqual(code, 1)
            self.assertIn("--final-eff: no such file", err)


class CootFilesTests(unittest.TestCase):
    def test_opens_a_model_and_its_maps_outside_a_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, maps, dictionary = root / "refined_001.pdb", root / "maps.mtz", root / "DU.cif"
            for path in (model, maps, dictionary):
                path.write_text("x\n")
            coot = root / "coot"
            coot.write_text("#!/bin/sh\n")
            calls = []
            view = launch_coot_files(
                model, maps, coot, working_directory=root / "run_001" / "CootGUI",
                source="AutoFDP f''-only round", dictionaries=[dictionary],
                launcher=lambda command, **kwargs: calls.append((command, kwargs)) or SimpleNamespace(pid=7),
            )
            command, kwargs = calls[0]
            self.assertEqual(command[1:3], ["--pdb", str(model.resolve())])
            self.assertIn("--auto", command)
            self.assertEqual(command[command.index("--auto") + 1], str(maps.resolve()))
            self.assertEqual(command[command.index("--dictionary") + 1], str(dictionary.resolve()))
            self.assertEqual(kwargs["cwd"], (root / "run_001" / "CootGUI").resolve())
            record = json.loads((view.working_directory / "launch.json").read_text())
            self.assertEqual((record["pid"], record["stage"]), (7, "files"))
            self.assertIn("ANOM", (view.working_directory / "load_anomalous_map.py").read_text())
            with self.assertRaisesRegex(CootViewError, "map coefficients is missing"):
                launch_coot_files(model, root / "none.mtz", coot,
                                  working_directory=root / "w", source="x")



if __name__ == "__main__":
    unittest.main()
