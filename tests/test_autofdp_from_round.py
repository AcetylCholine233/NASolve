"""``nasolve autofdp ROUND``: AutoFDP from an AutoRefine round folder.

The round's ``refined_NNN.pdb`` is the starting model, the reflection file its
``refined_NNN.eff`` refined against is the MTZ, and that .eff is the final
round's (user's request, 2026-11-19). phenix.refine, phenix.ready_set,
phenix.mtz.dump, phenix.python and Coot are simulated.
"""

import json
import tempfile
import unittest
from pathlib import Path

from nasolve.autofdp import AutoFDPError, read_autorefine_round

from .test_autofdp_final import FinalCommandTests, make_final_refine
from .test_autofdp_round import _plan, make_inputs


def write_round(model: Path, mtz: Path, *, labels: str = "F(+),SIGF(+),F(-),SIGF(-)",
                array_type: str = "amplitude", extra_array: str = "",
                report: dict | None = None) -> Path:
    """What AutoRefine leaves in a round: phenix.refine's refined_001.eff
    (output format: label groups with array_type, the round's own
    autorefine.params among phil_files, PostMR's model as the input model) and
    report.json."""
    round_dir = model.parent
    run = round_dir.parents[1]
    restraints = run / "PostMR" / "Restraints"
    restraints.mkdir(parents=True, exist_ok=True)
    (restraints / "DU.cif").write_text("data_comp_DU\n")
    (restraints / "pairs.phil").write_text("refinement { }\n")
    postmr_model = run / "PostMR" / "Model" / "readyset_model.pdb"
    postmr_model.parent.mkdir(parents=True, exist_ok=True)
    postmr_model.write_text(model.read_text())
    params = round_dir / "autorefine.params"
    params.write_text(
        "data_manager {\n  miller_array {\n"
        f'    file = "{mtz}"\n'
        f'    labels {{\n      name = "{labels}"\n    }}\n'
        '    labels {\n      name = "FreeR_flag"\n    }\n  }\n}\n'
        "refinement {\n  refine {\n    anomalous_scatterers {\n      group {\n"
        '        selection = "chain B and resid 101 and name AG"\n'
        "        refine = *f_prime *f_double_prime\n      }\n    }\n  }\n}\n"
    )
    eff = round_dir / "refined_001.eff"
    eff.write_text(
        "data_manager {\n"
        f'  restraint_files = "{restraints / "DU.cif"}"\n'
        f'  default_restraint = "{restraints / "DU.cif"}"\n'
        f'  phil_files = "{restraints / "pairs.phil"}"\n'
        f'  phil_files = "{params}"\n'
        f'  default_phil = "{restraints / "pairs.phil"}"\n'
        "  model {\n"
        f'    file = "{postmr_model}"\n'
        "    type = *x_ray neutron electron reference\n"
        "  }\n"
        f'  default_model = "{postmr_model}"\n'
        "  miller_array {\n"
        f'    file = "{mtz}"\n'
        "    labels {\n"
        f'      name = "{labels}"\n'
        "      type = *x_ray neutron electron\n"
        f"      array_type = unknown {'*' if array_type == 'amplitude' else ''}amplitude bool complex "
        "hendrickson_lattman \\\n"
        f"                   integer {'*' if array_type == 'intensity' else ''}intensity nonsense\n"
        "    }\n"
        "    labels {\n"
        '      name = "FreeR_flag"\n'
        "      type = *x_ray neutron electron\n"
        "      array_type = unknown amplitude bool complex hendrickson_lattman \\\n"
        "                   *integer intensity nonsense\n"
        "    }\n"
        f'    user_selected_labels = "{labels}"\n'
        '    user_selected_labels = "FreeR_flag"\n'
        "  }\n"
        + extra_array +
        "}\n"
        "refinement {\n"
        "  refine {\n"
        "    strategy = *individual_sites *individual_sites_real_space rigid_body \\\n"
        "               individual_adp *group_adp tls *occupancies *group_anomalous den\n"
        "    anomalous_scatterers {\n"
        "      group {\n"
        '        selection = "chain B and resid 101 and name AG"\n'
        "        f_prime = 0\n"
        "        f_double_prime = 0\n"
        "        refine = *f_prime *f_double_prime\n"
        "      }\n"
        "    }\n"
        "  }\n"
        "}\n"
        "output {\n"
        '  prefix = "refined"\n'
        "  serial = 1\n"
        "}\n"
    )
    (round_dir / "report.json").write_text(json.dumps(report if report is not None else {
        "checkpoint": "refine-002", "status": "AUTOREFINE_READY",
        "outputs": {"model": str(model)},
    }))
    return round_dir


class ReadRoundTests(unittest.TestCase):
    def test_model_mtz_and_eff_come_from_the_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            round_dir = write_round(model, mtz)
            source = read_autorefine_round(round_dir)
            self.assertEqual(source.directory, round_dir.resolve())
            self.assertEqual(source.model, model.resolve())
            self.assertEqual(source.eff, round_dir.resolve() / "refined_001.eff")
            self.assertEqual(source.mtz, mtz.resolve())
            # Free-R flags are not data.
            self.assertEqual(source.observation_labels, ("F(+),SIGF(+),F(-),SIGF(-)",))
            self.assertEqual((source.checkpoint, source.status), ("refine-002", "AUTOREFINE_READY"))

            plan = _plan(root, source.model, source.mtz, final_eff=source.eff, source=source)
            self.assertIs(plan.source, source)
            self.assertEqual(plan.warnings, ())
            self.assertEqual(plan.final_eff.source, source.eff)
            self.assertEqual(plan.final_eff.warnings, ())   # the same MTZ
            self.assertEqual(plan.final_eff.job_phil_files, (str(round_dir / "autorefine.params"),))
            with self.assertRaisesRegex(AutoFDPError, "must be the AutoRefine round's"):
                _plan(root, source.model, source.mtz, source=source)

    def test_a_separate_free_r_file_or_phases_are_not_the_data(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            other = mtz.parent / "phases.mtz"
            other.write_bytes(b"phases")
            round_dir = write_round(model, mtz, extra_array=(
                "  miller_array {\n"
                f'    file = "{other}"\n'
                "    labels {\n"
                '      name = "HLAM,HLBM,HLCM,HLDM"\n'
                "      array_type = unknown amplitude bool complex *hendrickson_lattman \\\n"
                "                   integer intensity nonsense\n"
                "    }\n"
                "  }\n"
            ))
            self.assertEqual(read_autorefine_round(round_dir).mtz, mtz.resolve())

    def test_intensities_in_the_eff_are_reported_against_autofdps_choice(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            source = read_autorefine_round(write_round(
                model, mtz, labels="I(+),SIGI(+),I(-),SIGI(-)", array_type="intensity"))
            self.assertEqual(source.observation_labels, ("I(+),SIGI(+),I(-),SIGI(-)",))
            plan = _plan(root, source.model, source.mtz, final_eff=source.eff, source=source)
            self.assertEqual(plan.warnings, (
                "refined_001.eff refined against I(+),SIGI(+),I(-),SIGI(-); the f'' rounds use "
                "F(+),SIGF(+),F(-),SIGF(-) (the first anomalous array in the MTZ), the final "
                "round the .eff's",
            ))

    def test_what_is_not_a_usable_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            round_dir = write_round(model, mtz)
            run = round_dir.parents[1]
            with self.assertRaisesRegex(AutoFDPError, "is not an AutoRefine round; give one of its rounds: "
                                        + str(round_dir.resolve()).replace("\\", "\\\\")):
                read_autorefine_round(run)
            with self.assertRaisesRegex(AutoFDPError, "is not a folder"):
                read_autorefine_round(model)
            empty = run / "AutoRefine" / "round_003"
            empty.mkdir()
            with self.assertRaisesRegex(AutoFDPError, "No refined_NNN.pdb in .*did not finish"):
                read_autorefine_round(empty)

            eff_text = (round_dir / "refined_001.eff").read_text()
            cases = {
                "No refined_001.eff beside refined_001.pdb": None,
                "names more than one reflection file with X-ray data": eff_text.replace(
                    "}\nrefinement {",
                    f'  miller_array {{\n    file = "{root / "other.mtz"}"\n'
                    '    labels {\n      name = "F,SIGF"\n'
                    "      array_type = unknown *amplitude intensity\n    }\n  }\n}\nrefinement {", 1),
                "names no reflection file with X-ray data": eff_text.replace("*amplitude", "amplitude"),
                "The reflection file refined_001.eff refined against does not exist": eff_text.replace(
                    str(mtz), str(root / "moved.mtz")),
            }
            for message, text in cases.items():
                with self.subTest(message=message):
                    if text is None:
                        (round_dir / "refined_001.eff").unlink()
                    else:
                        self.assertNotEqual(text, eff_text)
                        (round_dir / "refined_001.eff").write_text(text)
                    with self.assertRaisesRegex(AutoFDPError, message):
                        read_autorefine_round(round_dir)
            (round_dir / "refined_001.eff").write_text(eff_text)
            (round_dir / "report.json").write_text(json.dumps({"outputs": {"model": "refined_009.pdb"}}))
            with self.assertRaisesRegex(AutoFDPError, "records the model refined_009.pdb, not refined_001.pdb"):
                read_autorefine_round(round_dir)
            (round_dir / "refined_002.pdb").write_text(model.read_text())
            with self.assertRaisesRegex(AutoFDPError, "more than one refined model"):
                read_autorefine_round(round_dir)


class RoundCommandTests(unittest.TestCase):
    _main = FinalCommandTests._main

    def test_autofdp_runs_from_the_round_and_ends_with_its_eff(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            round_dir = write_round(model, mtz)
            code, out, err, launched = self._main(root, [str(round_dir)], make_final_refine(root))
            self.assertEqual(code, 0, err)
            eff = round_dir.resolve() / "refined_001.eff"
            for text in (
                f"AutoRefine round: {round_dir.resolve()} (refine-002, AUTOREFINE_READY)",
                f"Model: {model.resolve()}",
                f"MTZ: {mtz.resolve()}",
                f"Final round settings (the AutoRefine round's .eff): {eff}",
                f"left out: {round_dir / 'autorefine.params'} (a refinement job's own parameters",
                "AutoFDP final round and B confirmation (final_round, settings of refined_001.eff; ",
                "AUTOFDP_CONVERGED",
            ):
                self.assertIn(text, out)
            run = mtz.parent.resolve() / "AutoFDP" / "run_001"
            report = json.loads((run / "report.json").read_text())
            recorded = report["inputs"]["autorefine_round"]
            self.assertEqual(recorded["directory"], {
                "anchor": "dataset", "relative_path": "AutoMR/run_001/AutoRefine/round_002"})
            self.assertEqual((recorded["checkpoint"], recorded["status"]), ("refine-002", "AUTOREFINE_READY"))
            self.assertEqual(recorded["eff"]["relative_path"], "AutoMR/run_001/AutoRefine/round_002/refined_001.eff")
            self.assertEqual(recorded["eff_observation_labels"], ["F(+),SIGF(+),F(-),SIGF(-)"])
            self.assertEqual(report["inputs"]["model"]["relative_path"],
                             "AutoMR/run_001/AutoRefine/round_002/refined_001.pdb")
            self.assertEqual(report["inputs"]["mtz"]["relative_path"], "staraniso_alldata-unique.mtz")
            self.assertEqual(report["settings"]["final_eff"], str(eff))
            self.assertEqual(report["warnings"], [])
            self.assertIn(f"AutoRefine round: {round_dir.resolve()} (refine-002)",
                          (run / "autofdp.log").read_text())
            launched.assert_called_once()

            # The model-and-MTZ form is unchanged and records no round.
            code, out, err, _ = self._main(root, [str(mtz), str(model), "--no-coot"], make_final_refine(root))
            self.assertEqual(code, 0, err)
            self.assertNotIn("AutoRefine round:", out)
            report = json.loads((mtz.parent / "AutoFDP" / "run_002" / "report.json").read_text())
            self.assertIsNone(report["inputs"]["autorefine_round"])

    def test_round_problems_on_the_command_line(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            round_dir = write_round(model, mtz)
            other = root / "manual.eff"
            other.write_text("refinement { }\n")
            code, _, err, _ = self._main(root, [str(round_dir), "--final-eff", str(other)],
                                         make_final_refine(root))
            self.assertEqual(code, 1)
            self.assertIn("--final-eff cannot be given with an AutoRefine round", err)
            code, _, err, _ = self._main(root, [str(round_dir.parents[1])], make_final_refine(root))
            self.assertEqual(code, 1)
            self.assertIn("is not an AutoRefine round; give one of its rounds", err)
            code, _, err, _ = self._main(root, [str(model), str(mtz), str(round_dir)], make_final_refine(root))
            self.assertEqual(code, 1)
            self.assertIn("Give one model (.pdb) and one reflection file (.mtz)", err)
            self.assertFalse((mtz.parent / "AutoFDP").exists())


if __name__ == "__main__":
    unittest.main()
