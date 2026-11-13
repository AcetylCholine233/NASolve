"""AutoRefine and RefineDoctor take the X-ray wavelength from the MTZ header.

Regression: the wavelength used for calculated (Henke) f'/f'' came only from
an AutoSol record, so a metal structure refined without AutoSol (7CU_Ag,
Phenix 2.2.1-6174) reported no calculated f'' and RefineDoctor skipped every
fixed/fdp-only scattering recipe. Simulated with fake Phenix executables.
"""

import json
import tempfile
import unittest
from pathlib import Path

from nasolve.autorefine import (
    AutoRefineError,
    build_reflection_plan,
    calculated_anomalous_groups,
    execute_autorefine,
)
from nasolve.refine_doctor import execute_refine_doctor

from .test_autorefine import PHENIX_21, make_mtz_dump, make_refine, make_refine_run

MTZ = "MTZ header (phenix.mtz.dump)"


def _plan(run: Path, dump: Path):
    report = json.loads((run / "report.json").read_text())
    observations = Path(report["inputs"]["reflections"])
    return report, build_reflection_plan(report, observations, dump)


class ReflectionPlanWavelengthTests(unittest.TestCase):
    def test_plan_reads_the_mtz_header_without_autosol(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            _, plan = _plan(run, make_mtz_dump(root, wavelength=1.00743))
            self.assertEqual(plan.wavelength, 1.00743)
            self.assertEqual(plan.wavelength_source["method"], MTZ)
            self.assertEqual(plan.wavelength_source["mtz_dataset"], "1")
            self.assertIsNone(plan.wavelength_source["autosol_wavelength"])

    def test_mtz_header_wins_over_a_disagreeing_autosol_record(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)  # AutoSol recorded 1.377618 A
            _, plan = _plan(run, make_mtz_dump(root, wavelength=1.37762))
            self.assertEqual(plan.wavelength, 1.37762)
            self.assertEqual(plan.wavelength_source["warnings"], [])
            _, plan = _plan(run, make_mtz_dump(root, wavelength=0.9795))
            self.assertEqual(plan.wavelength, 0.9795)
            self.assertEqual(plan.wavelength_source["method"], MTZ)
            self.assertIn("MTZ header value is used", plan.wavelength_source["warnings"][0])

    def test_autosol_record_is_only_a_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            _, plan = _plan(run, make_mtz_dump(root))  # header has no dataset blocks
            self.assertEqual(plan.wavelength, 1.377618)
            self.assertEqual(plan.wavelength_source["method"], "autosol.wavelength")
            self.assertIn("wavelength", plan.wavelength_source["mtz_error"])

            bare = make_refine_run(root / "bare", autosol=False)
            report, plan = _plan(bare, make_mtz_dump(root / "bare"))
            self.assertIsNone(plan.wavelength)
            self.assertIsNone(plan.wavelength_source["method"])
            with self.assertRaisesRegex(AutoRefineError, "MTZ header has none"):
                calculated_anomalous_groups(report, make_refine(root / "bare"), None, plan)

    def test_calculated_groups_record_the_mtz_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            report, plan = _plan(run, make_mtz_dump(root, wavelength=1.37762))
            groups = calculated_anomalous_groups(report, make_refine(root), None, plan)
            group = groups["chain B and resid 4 and name I"]
            self.assertEqual(group["wavelength"], 1.37762)
            self.assertEqual(group["wavelength_source"], MTZ)
            self.assertEqual(group["f_double_prime"], 5.75)


class RefinementWavelengthTests(unittest.TestCase):
    def test_round_without_autosol_reports_calculated_f_double_prime(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            result = execute_autorefine(
                run, make_refine(root), make_mtz_dump(root, wavelength=1.37762),
                phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
                processor_count=8,
            )
            site = result.statistics["anomalous_scatterers"][0]
            self.assertEqual(site["wavelength"], 1.37762)
            self.assertEqual(site["wavelength_source"], MTZ)
            self.assertAlmostEqual(site["calculated_f_double_prime"], 5.750252723693848)
            self.assertIsNotNone(site["apparent_anomalous_occupancy"])
            payload = json.loads(result.report_path.read_text())
            self.assertEqual(payload["inputs"]["wavelength"], 1.37762)
            self.assertEqual(payload["inputs"]["wavelength_source"]["method"], MTZ)

    def test_doctor_runs_explicit_scattering_recipes_without_autosol(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root, autosol=False)
            refine = make_refine(root, final_work=.2, final_free=.19)
            result = execute_refine_doctor(
                run, refine, make_mtz_dump(root, wavelength=1.37762),
                phenix_version=PHENIX_21, max_trials=3,
            )
            explicit = []
            for trial in result.trials:
                payload = json.loads(trial.report_path.read_text())
                if payload["refinement"]["anomalous_parameter_mode"] in {"fixed", "fdp-only"}:
                    explicit.append(payload["refinement"]["anomalous_groups"])
            self.assertTrue(explicit, "no fixed/fdp-only scattering trial ran")
            for groups in explicit:
                group = groups["chain B and resid 4 and name I"]
                self.assertEqual(group["wavelength"], 1.37762)
                self.assertEqual(group["wavelength_source"], MTZ)


if __name__ == "__main__":
    unittest.main()
