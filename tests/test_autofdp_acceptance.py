"""When a refined f'' matches theory (user's request, 2026-11-23).

The default is +/-10% of the theoretical f'' of the site's element, so a weak
anomalous scatterer such as Ag (f'' about 2 e- at 1 A) is judged on the same
relative scale as Hg (about 10 e-). It was a fixed +/-1 e- (``--tolerance``),
which for Ag accepted occupancies from about 0.5 to 1 and more.
``--acceptance_absolute E`` gives a fixed window of E electrons instead;
``--tolerance`` is gone.

phenix.refine, phenix.ready_set, phenix.mtz.dump and phenix.python are
simulated (theory: Ag f'' 2.03, Hg f'' 10.1).
"""

import json
import tempfile
import unittest
from pathlib import Path

from nasolve.autofdp import (
    ACCEPTANCE_FRACTION,
    STATUS_CONVERGED,
    STATUS_REVIEW,
    AutoFDPError,
    acceptance_window,
)

from . import test_autofdp_round as round_tests
from .test_autofdp_round import MERCURY, SILVER, _execute, _plan, make_fdp_refine, make_inputs


class AcceptanceWindowTests(unittest.TestCase):
    def test_ten_percent_of_theory_unless_absolute(self):
        self.assertEqual(ACCEPTANCE_FRACTION, 0.10)
        self.assertAlmostEqual(acceptance_window(2.03), 0.203)
        self.assertAlmostEqual(acceptance_window(10.1), 1.01)
        self.assertEqual(acceptance_window(2.03, absolute=1.0), 1.0)
        for bad in (0, -1, float("nan"), float("inf")):
            with self.subTest(absolute=bad), self.assertRaisesRegex(
                    AutoFDPError, "--acceptance_absolute must be a positive number of electrons"):
                acceptance_window(2.03, absolute=bad)


class DefaultAcceptanceTests(unittest.TestCase):
    def test_ag_off_by_more_than_ten_percent_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)              # Ag B:101 at occupancy 0.90
            refine = make_fdp_refine(root, {SILVER: 2.6})   # +0.57, within the old 1 e-
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine,
                              max_iterations=0)
            site = result.round.sites[0]
            self.assertAlmostEqual(site.acceptance, 0.203)
            self.assertFalse(site.within_acceptance)
            self.assertEqual(result.status, STATUS_REVIEW)
            self.assertEqual(result.round.message, "1 of 1 site(s) outside +/-10% of theory")

            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["schema_version"], 3)
            self.assertEqual(report["settings"]["acceptance"],
                             {"rule": "relative", "fraction": 0.1, "absolute_electrons": None})
            self.assertNotIn("tolerance", report["settings"])
            entry = report["rounds"][0]["sites"][0]
            self.assertAlmostEqual(entry["acceptance"], 0.203)
            self.assertFalse(entry["within_acceptance"])
            self.assertNotIn("within_tolerance", entry)
            log = result.log_path.read_text()
            self.assertIn("Acceptance: refined f'' within +/-10% of theory (AG +/-0.20 e-)", log)

    def test_each_element_gets_its_own_window(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)
            # Ag +0.17 (window 0.203) and Hg -0.90 (window 1.01): both accepted.
            refine = make_fdp_refine(root, {SILVER: 2.2, MERCURY: 9.2})
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine)
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)
            silver, mercury = result.round.sites
            self.assertAlmostEqual(silver.acceptance, 0.203)
            self.assertAlmostEqual(mercury.acceptance, 1.01)
            self.assertTrue(silver.within_acceptance and mercury.within_acceptance)
            self.assertEqual(result.outcomes[0].reason,
                             "within +/-10% (0.20 e-) of theory without adjustment")
            self.assertEqual(result.message, "every site within +/-10% of theory after "
                                             "0 occupancy adjustment(s) in 1 f'' round(s)")

    def test_acceptance_absolute_gives_a_fixed_window(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_fdp_refine(root, {SILVER: 2.6})
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine,
                              acceptance_absolute=1.0)
            site = result.round.sites[0]
            self.assertEqual((site.acceptance, site.within_acceptance), (1.0, True))
            self.assertEqual(result.status, STATUS_CONVERGED)
            self.assertEqual(result.outcomes[0].reason,
                             "within +/-1 e- of theory without adjustment")
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["settings"]["acceptance"],
                             {"rule": "absolute", "fraction": None, "absolute_electrons": 1.0})
            self.assertIn("Acceptance: refined f'' within +/-1 e- of theory",
                          result.log_path.read_text())


class AcceptanceCommandTests(unittest.TestCase):
    _main = round_tests.AutoFDPCommandTests._main

    def test_default_and_acceptance_absolute(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            code, out, err, _ = self._main(root, [str(model), str(mtz), "--plan"])
            self.assertEqual(code, 0, err)
            self.assertIn("per site outside +/-10% of theory (AG +/-0.20 e-)", out)

            refine = make_fdp_refine(root, {SILVER: 2.6})
            for spelling in ("--acceptance_absolute", "--acceptance-absolute"):
                with self.subTest(spelling=spelling):
                    code, out, err, _ = self._main(root, [str(model), str(mtz), spelling, "1"],
                                                   refine=refine)
                    self.assertEqual(code, 0, err)
                    row = next(line for line in out.splitlines() if line.startswith("B:101 AG"))
                    self.assertIn("+0.57  within +/-1.00", row)
                    self.assertIn("every site already within +/-1 e- of theory", out)

            # Default: +0.57 is outside +/-0.20 e- for Ag; review after round_001.
            code, out, err, _ = self._main(root, [str(model), str(mtz), "--max-iterations", "0"],
                                           refine=refine)
            self.assertEqual(code, 2, err)
            row = next(line for line in out.splitlines() if line.startswith("B:101 AG"))
            self.assertIn("+0.57  outside +/-0.20", row)

    def test_tolerance_is_gone_and_bad_values_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            with self.assertRaises(SystemExit) as stopped:
                self._main(root, [str(model), str(mtz), "--tolerance", "0.2"])
            self.assertEqual(stopped.exception.code, 2)
            code, _, err, _ = self._main(root, [str(model), str(mtz), "--acceptance_absolute", "0"])
            self.assertEqual(code, 1)
            self.assertIn("--acceptance_absolute must be a positive number of electrons", err)
            self.assertFalse((mtz.parent / "AutoFDP").exists())


if __name__ == "__main__":
    unittest.main()
