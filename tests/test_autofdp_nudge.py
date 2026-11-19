"""``nasolve autofdp --nudge Q`` (user's request, 2026-11-19): the heavy atom
starts at the occupancy AutoRefine left it at plus Q (first written as minus
Q; the user corrected it the same day). For now one metal site
is nudged: the model's only one, or the one --site names.

phenix.refine, phenix.ready_set, phenix.mtz.dump and phenix.python are
simulated.
"""

import json
import tempfile
import unittest
from pathlib import Path

from nasolve.autofdp import AutoFDPError, nudge_targets, prepare_start_model

from .test_autofdp_round import AutoFDPCommandTests, _execute, _plan, make_inputs


def occupancy(path: Path, atom: str = "AG") -> float:
    line = next(line for line in path.read_text().splitlines()
                if line.startswith("HETATM") and line[12:16].strip() == atom)
    return float(line[54:60])


class NudgeStartTests(unittest.TestCase):
    def test_start_occupancy_is_the_models_plus_the_nudge(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, _ = make_inputs(root)            # Ag B:101 at occupancy 0.90
            original = model.read_bytes()
            plan = prepare_start_model(model, nudge=0.05)
            site = plan.sites[0]
            self.assertEqual((site.occupancy, site.occupancy_start, site.nudge), (0.9, 0.95, 0.05))
            written = prepare_start_model(model, nudge=0.05, output_directory=root / "run")
            self.assertEqual(occupancy(written.start_model), 0.95)
            self.assertEqual(written.source_copy.read_bytes(), original)   # never changed
            self.assertEqual(model.read_bytes(), original)
            # Only the occupancy and B columns of the metal line change.
            start = written.start_model.read_text().splitlines()
            source = written.source_copy.read_text().splitlines()
            changed = [(a, b) for a, b in zip(source, start) if a != b]
            self.assertEqual(len(changed), 1)
            self.assertEqual(changed[0][0][:54], changed[0][1][:54])
            self.assertEqual(changed[0][0][66:], changed[0][1][66:])

            # A negative nudge lowers it; without --nudge nothing changes.
            self.assertEqual(prepare_start_model(model, nudge=-0.1).sites[0].occupancy_start, 0.8)
            plain = prepare_start_model(model).sites[0]
            self.assertEqual((plain.occupancy_start, plain.nudge), (0.9, None))

    def test_nudge_limits(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, _ = make_inputs(root)
            for nudge, floor, message in (
                (0.2, 0.05, "above 1"),
                (-0.88, 0.05, "below the occupancy floor 0.05"),
                (-0.5, 0.5, "below the occupancy floor 0.5"),
                (float("nan"), 0.05, "--nudge must be a number"),
            ):
                with self.subTest(nudge=nudge, floor=floor):
                    with self.assertRaisesRegex(AutoFDPError, message):
                        prepare_start_model(model, nudge=nudge, occupancy_floor=floor)

    def test_one_site_for_now(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, _ = make_inputs(root, mercury=True)   # Ag B:101 and Hg C:1
            with self.assertRaisesRegex(
                AutoFDPError, r"--nudge sets the starting occupancy of one metal site for now; "
                              r"this model has 2 \(B:101, C:1\): name it with --site"):
                prepare_start_model(model, nudge=0.1)
            sites = prepare_start_model(model, nudge=0.1, sites=["C:1"]).sites
            self.assertEqual([(site.site, site.occupancy, site.occupancy_start) for site in sites],
                             [("C:1", 0.5, 0.6)])
            # The selection is its own function, to be widened later.
            plan = prepare_start_model(model, sites=["B:101"])
            self.assertEqual([site.site for site in nudge_targets(plan.sites)], ["B:101"])


class NudgeRunTests(unittest.TestCase):
    def test_the_first_round_starts_from_the_nudged_occupancy_and_it_is_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            plan = _plan(root, model, mtz, nudge=0.05)
            self.assertEqual(plan.nudge, 0.05)
            self.assertEqual(plan.b_start.sites[0].occupancy_start, 0.95)
            result = _execute(plan, root, max_iterations=0)
            run = result.run_directory
            self.assertEqual(occupancy(run / "start_model.pdb"), 0.95)
            self.assertEqual(occupancy(run / "readyset_model.pdb"), 0.95)
            self.assertEqual(result.rounds[0].sites[0].occupancy, 0.95)
            report = json.loads((run / "report.json").read_text())
            self.assertEqual(report["settings"]["nudge"], 0.05)
            site = report["b_start"]["sites"][0]
            self.assertEqual((site["occupancy"], site["occupancy_start"], site["nudge"]), (0.9, 0.95, 0.05))
            self.assertIn("occupancy 0.90 -> 0.95 (--nudge 0.05)", (run / "autofdp.log").read_text())


class NudgeCommandTests(unittest.TestCase):
    _main = AutoFDPCommandTests._main

    def test_command_line(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            code, out, err, _ = self._main(root, [str(model), str(mtz), "--nudge", "0.05", "--plan"])
            self.assertEqual(code, 0, err)
            self.assertIn("B:101 AG (AG), occupancy 0.90 -> 0.95 (--nudge 0.05): B 379.99 -> 50.00", out)
            self.assertFalse((mtz.parent / "AutoFDP").exists())
            code, out, err, _ = self._main(root, [str(model), str(mtz), "--nudge", "-0.9"])
            self.assertEqual(code, 1)
            self.assertIn("below the occupancy floor 0.05", err)
            self.assertFalse((mtz.parent / "AutoFDP").exists())


if __name__ == "__main__":
    unittest.main()
