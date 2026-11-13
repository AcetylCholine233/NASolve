"""AutoFDP step 3: occupancy adjustment.

The stand-in phenix.refine answers the occupancy it is given: for each
anomalous group it reads the site's occupancy q from the input model and
reports f'' = a + S/q, which is how a real f''-only refinement behaves when
the anomalous signal fixes q x f''. So the loop, the proportional and secant
steps, the bounds and the per-site limit run end to end. phenix.refine,
phenix.ready_set, phenix.mtz.dump, phenix.python and Coot are simulated.
"""

import json
import stat
import tempfile
import unittest
from pathlib import Path

from nasolve.autofdp import (
    STATUS_CONVERGED,
    STATUS_FAILED,
    STATUS_REVIEW,
    STATUS_ROUND_DONE,
    AutoFDPError,
    propose_occupancy,
    write_site_occupancy,
    wrong_response,
)

from .helpers import pdb_record
from .test_autofdp_round import (
    MERCURY,
    SILVER,
    AutoFDPCommandTests,
    _execute,
    _plan,
    make_fdp_refine,
    make_inputs,
)


def make_scaling_refine(root: Path, signal: dict[str, object], *,
                        fail_below: float | None = None) -> Path:
    """Stand-in phenix.refine reporting f'' = a + S/q for each group in
    ``signal`` (selection -> (a, S)), with q read from the input model.

    A dict instead of (a, S) can bend the curve the way Hg h:1 of PdUHgPdU
    did: ``peak`` makes f'' fall again below that occupancy (a + S q/peak^2),
    ``trough`` makes it rise again above it; ``peak_while`` = (selection,
    occupancy) keeps the peak only while that other site is above the
    occupancy. ``fail_below``: exit 1 when a site's occupancy is below it."""
    make_fdp_refine(root)          # phenix.python with the theory: AG 2.03, HG 10.1
    executable = root / "phenix.refine"
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
        "lines = model.read_text().splitlines()\n"
        f"signal = {signal!r}\n"
        f"fail_below = {fail_below!r}\n"
        "def occupancy(selection):\n"
        "    chain, resid, name = re.fullmatch(\n"
        "        r'chain (\\S+) and resid (\\S+) and name (\\S+)', selection).groups()\n"
        "    for line in lines:\n"
        "        if (line.startswith(('ATOM', 'HETATM')) and line[21] == chain\n"
        "                and line[22:26].strip() == resid and line[12:16].strip() == name):\n"
        "            return float(line[54:60])\n"
        "    raise SystemExit('no atom for ' + selection)\n"
        "values = {selection: occupancy(selection) for selection in signal}\n"
        "if fail_below is not None and any(q < fail_below for q in values.values()):\n"
        "    print('Sorry: fixture failure')\n"
        "    raise SystemExit(1)\n"
        "shutil.copyfile(model, 'refined_001.pdb')\n"
        "Path('refined_001_map_coeffs.mtz').write_bytes(b'maps')\n"
        "print('Final: r_work = 0.2161 r_free = 0.2121 bonds = 0.012 angles = 1.82')\n"
        "for block in re.findall(r'group \\{(.*?)\\}', params, re.S):\n"
        "    selection = re.search(r'selection = \"(.*?)\"', block).group(1)\n"
        "    if selection not in signal:\n"
        "        continue\n"
        "    spec = signal[selection]\n"
        "    if not isinstance(spec, dict):\n"
        "        spec = {'a': spec[0], 'S': spec[1]}\n"
        "    offset, scale, q = spec.get('a', 0.0), spec['S'], values[selection]\n"
        "    peak, trough, condition = spec.get('peak'), spec.get('trough'), spec.get('peak_while')\n"
        "    if condition is not None and values[condition[0]] <= condition[1]:\n"
        "        peak = None\n"
        "    if peak is not None and q < peak:\n"
        "        value = offset + scale * q / peak ** 2\n"
        "    elif trough is not None and q > trough:\n"
        "        value = offset + scale * q / trough ** 2\n"
        "    else:\n"
        "        value = offset + scale / q\n"
        "    print('Anomalous scatterer group:')\n"
        "    print(f'  Selection: \"{selection}\"')\n"
        "    print('  Number of selected scatterers: 1')\n"
        "    print('  f_prime: ' + re.search(r'f_prime = (\\S+)', block).group(1))\n"
        "    print(f'  f_double_prime: {value:.3f}')\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def occupancy_of(path: Path, chain: str, resid: str, name: str) -> float:
    for line in path.read_text().splitlines():
        if (line.startswith(("ATOM", "HETATM")) and line[21] == chain
                and line[22:26].strip() == resid and line[12:16].strip() == name):
            return float(line[54:60])
    raise AssertionError(f"no {chain}:{resid} {name} in {path}")


class ProposeOccupancyTests(unittest.TestCase):
    def test_first_step_is_proportional_in_the_direction_of_theory(self):
        # f'' below theory: the site holds less metal than modelled -> lower q.
        step = propose_occupancy([(0.9, 1.167)], 2.03)
        self.assertEqual((step.occupancy, step.method, step.bound), (0.52, "proportional", None))
        # f'' above theory -> raise q.
        self.assertEqual(propose_occupancy([(0.5, 3.0)], 2.03).occupancy, 0.74)

    def test_later_steps_use_a_secant_on_inverse_occupancy(self):
        # f'' = 0.5 + 0.6/q: the proportional step from 0.90 lands at 0.52;
        # the secant through both rounds solves 0.5 + 0.6/q = 2.03 exactly.
        step = propose_occupancy([(0.9, 1.167), (0.52, 1.654)], 2.03)
        self.assertEqual((step.occupancy, step.method), (0.39, "secant"))
        # Rounds at an unchanged occupancy (another site was adjusted) are
        # skipped: the secant uses the last round with another occupancy.
        step = propose_occupancy([(0.9, 1.1), (0.52, 1.6), (0.52, 1.62)], 2.03)
        self.assertEqual((step.occupancy, step.method), (0.39, "secant"))

    def test_secant_falls_back_to_proportional(self):
        # f'' fell as occupancy fell: no usable secant.
        step = propose_occupancy([(0.9, 1.0), (0.5, 0.9)], 2.03)
        self.assertEqual((step.occupancy, step.method), (0.22, "proportional"))
        self.assertIn("secant not used", step.note)

    def test_step_size_limits(self):
        step = propose_occupancy([(0.8, -0.3)], 2.03)   # no positive f'': at most 4-fold down
        self.assertEqual(step.occupancy, 0.2)
        self.assertEqual(step.note, "step limited to 4-fold")
        # 0.5 x 2.025 / 2.03 rounds back to 0.50: at least one 0.01 step.
        self.assertEqual(propose_occupancy([(0.5, 2.025)], 2.03).occupancy, 0.49)

    def test_bounds(self):
        step = propose_occupancy([(0.9, 3.0)], 2.03)
        self.assertEqual((step.occupancy, step.bound, step.note), (1.0, "ceiling", "stopped at full occupancy"))
        step = propose_occupancy([(1.0, 3.0)], 2.03)
        self.assertEqual((step.occupancy, step.bound, step.note),
                         (None, "ceiling", "f'' above theory at full occupancy"))
        step = propose_occupancy([(0.06, 0.5)], 2.03)
        self.assertEqual((step.occupancy, step.bound), (0.05, "floor"))
        step = propose_occupancy([(0.05, 0.5)], 2.03)
        self.assertEqual((step.occupancy, step.bound, step.note),
                         (None, "floor", "f'' below theory at the occupancy floor"))
        self.assertEqual(propose_occupancy([(0.9, 0.5)], 2.03, floor=0.3).occupancy, 0.3)

    def test_unusable_input(self):
        for history, theory, message in (
            ([], 2.03, "No f'' round"),
            ([(0.9, None)], 2.03, "no refined f''"),
            ([(0.9, 1.0)], 0.0, "must be a positive number"),
        ):
            with self.subTest(history=history, theory=theory):
                with self.assertRaisesRegex(AutoFDPError, message):
                    propose_occupancy(history, theory)
        with self.assertRaisesRegex(AutoFDPError, "occupancy floor"):
            propose_occupancy([(0.9, 1.0)], 2.03, floor=1.0)


class WrongResponseTests(unittest.TestCase):
    def test_fdp_falling_with_occupancy_returns_the_occupancy_with_the_highest_fdp(self):
        # PdUHgPdU AutoFDP/run_003, Hg h:1 (live, Phenix 2.2.1): f'' rose from
        # 1.00 to 0.34, then fell from 0.34 to 0.20.
        self.assertIsNone(wrong_response([(1.0, 3.44426), (0.34, 6.6448)]))
        response = wrong_response([(1.0, 3.44426), (0.34, 6.6448), (0.2, 6.38852)])
        self.assertEqual(
            (response.lowered, response.occupancy_before, response.occupancy_after,
             response.best_occupancy, response.best_f_double_prime),
            (True, 0.34, 0.2, 0.34, 6.6448),
        )

    def test_mirror_case_and_rounds_without_an_occupancy_change(self):
        response = wrong_response([(0.5, 12.12), (0.6, 13.447)])   # raised, f'' rose
        self.assertEqual((response.lowered, response.best_occupancy), (False, 0.5))
        self.assertIsNone(wrong_response([(0.5, 12.12), (0.6, 11.0)]))
        # Another site was adjusted: same occupancy, so no verdict.
        self.assertIsNone(wrong_response([(0.3, 7.4), (0.3, 7.0)]))
        self.assertIsNone(wrong_response([(0.5, 4.0)]))


class WriteSiteOccupancyTests(unittest.TestCase):
    def test_only_the_site_occupancy_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "in.pdb"
            source.write_text(
                pdb_record("ATOM", 1, "P", "DC", "A", 12, x=1.0, y=2.0, z=3.0)
                + pdb_record("HETATM", 2, "AG", "AG", "B", 101, element="AG",
                             x=4.0, y=5.0, z=6.0, occupancy=0.9, b_factor=50.0)
                + "END\n"
            )
            destination = root / "out.pdb"
            write_site_occupancy(source, destination, "B:101", "AG", None, 0.52)
            before, after = source.read_text().splitlines(), destination.read_text().splitlines()
            self.assertEqual(before[0], after[0])
            self.assertEqual(after[1][54:60], "  0.52")
            self.assertEqual(before[1][:54] + before[1][60:], after[1][:54] + after[1][60:])
            self.assertEqual(before[2:], after[2:])
            with self.assertRaisesRegex(AutoFDPError, "matches 0 atoms"):
                write_site_occupancy(source, destination, "B:102", "AG", None, 0.5)


class OccupancyLoopTests(unittest.TestCase):
    def test_site_converges_by_proportional_then_secant_step(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)            # Ag B:101 at occupancy 0.90
            refine = make_scaling_refine(root, {SILVER: (0.5, 0.6)})
            result = _execute(_plan(root, model, mtz, refine=refine), root,
                              refine=refine, tolerance=0.1)
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)
            run = result.run_directory
            self.assertEqual([item.directory.name for item in result.rounds],
                             ["round_001", "round_002", "round_003"])
            self.assertEqual(
                [(a.round, a.occupancy_before, a.occupancy_after, a.method) for a in result.adjustments],
                [(2, 0.9, 0.52, "proportional"), (3, 0.52, 0.39, "secant")],
            )
            self.assertEqual([item.sites[0].f_double_prime_refined for item in result.rounds],
                             [1.167, 1.654, 2.038])
            self.assertTrue(result.round.sites[0].within_tolerance)
            self.assertIn("after 2 occupancy adjustment(s) in 3 f'' round(s)", result.message)

            # Each round gets its own copy of the previous model with only the
            # site's occupancy changed, and the same parameters.
            for number, occupancy in ((2, 0.52), (3, 0.39)):
                received = json.loads((run / f"round_00{number}" / "received_args.json").read_text())
                self.assertIn(str(run / f"round_00{number}" / "input_model.pdb"), received)
                self.assertEqual(occupancy_of(run / f"round_00{number}" / "input_model.pdb",
                                              "B", "101", "AG"), occupancy)
            first = (run / "readyset_model.pdb").read_text().splitlines()
            final = result.final_model.read_text().splitlines()
            self.assertEqual(result.final_model, run / "final_model.pdb")
            changed = [index for index, (a, b) in enumerate(zip(first, final)) if a != b]
            self.assertEqual(len(changed), 1)
            self.assertEqual(final[changed[0]][54:66], "  0.39 50.00")   # B from step 1 kept
            params = [(run / f"round_00{n}" / "autofdp.params").read_text() for n in (1, 2, 3)]
            self.assertEqual(params[0], params[1])
            self.assertEqual(params[0], params[2])
            self.assertIn("strategy = *group_anomalous\n", params[0])
            self.assertIn("f_double_prime = 0.0", params[0])

            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["status"], STATUS_CONVERGED)
            self.assertEqual(report["settings"]["max_iterations"], 5)
            self.assertEqual(report["steps"]["occupancy_adjustment"], "done")
            self.assertEqual(report["steps"]["fdp_rounds"], "3 run")
            entry = report["rounds"][2]
            self.assertEqual(
                (entry["adjusted_site"], entry["occupancy_before"], entry["occupancy_after"],
                 entry["method"], entry["input_model"]["relative_path"]),
                (SILVER, 0.52, 0.39, "secant", "round_003/input_model.pdb"),
            )
            self.assertEqual(entry["sites"][0]["occupancy"], 0.39)
            site = report["occupancy"]["sites"][0]
            self.assertEqual((site["occupancy_start"], site["occupancy_final"], site["adjustments"],
                              site["converged"]), (0.9, 0.39, 2, True))
            self.assertEqual(site["reason"], "within +/-0.1 e- of theory after 2 adjustment(s)")
            self.assertEqual(
                [(item["round"], item["occupancy"], item["f_double_prime"]) for item in site["history"]],
                [(1, 0.9, 1.167), (2, 0.52, 1.654), (3, 0.39, 2.038)],
            )
            self.assertEqual(report["final_model"]["relative_path"], "final_model.pdb")
            log = result.log_path.read_text()
            self.assertIn("Step 3 (occupancy adjustment", log)
            self.assertIn("round_003: B:101 AG occupancy 0.52 -> 0.39 (secant)", log)

    def test_fdp_above_theory_at_full_occupancy_goes_to_review(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_scaling_refine(root, {SILVER: (0.0, 2.436)})   # needs q = 1.2
            result = _execute(_plan(root, model, mtz, refine=refine), root,
                              refine=refine, tolerance=0.1)
            self.assertEqual(result.status, STATUS_REVIEW)
            self.assertEqual([(a.occupancy_after, a.bound) for a in result.adjustments],
                             [(1.0, "ceiling")])
            self.assertEqual(result.outcomes[0].reason, "f'' above theory at full occupancy")
            self.assertEqual(result.message, "B:101 AG: f'' above theory at full occupancy")
            self.assertEqual(occupancy_of(result.final_model, "B", "101", "AG"), 1.0)
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["steps"]["occupancy_adjustment"], "review")

    def test_fdp_below_theory_at_the_occupancy_floor_goes_to_review(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_scaling_refine(root, {SILVER: (0.0, 0.0406)})  # needs q = 0.02
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine)
            self.assertEqual(result.status, STATUS_REVIEW)
            self.assertEqual([a.occupancy_after for a in result.adjustments], [0.23, 0.06, 0.05])
            self.assertEqual(result.outcomes[0].reason, "f'' below theory at the occupancy floor")
            raised = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine,
                              occupancy_floor=0.3)
            self.assertEqual([a.occupancy_after for a in raised.adjustments], [0.3])
            self.assertEqual(raised.outcomes[0].occupancy_final, 0.3)

    def test_each_site_gets_at_most_max_iterations_adjustments(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_scaling_refine(root, {SILVER: (0.5, 0.6)})
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine,
                              tolerance=0.1, max_iterations=1)
            self.assertEqual(result.status, STATUS_REVIEW)
            self.assertEqual(len(result.rounds), 2)
            self.assertEqual(
                result.outcomes[0].reason,
                "still -0.38 e- from theory after 1 adjustment(s) (--max-iterations)",
            )

    def test_site_furthest_from_theory_is_adjusted_first(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)   # Ag 0.90 and Hg 0.50
            # Ag needs q = 0.60 (f'' 1.35, -0.68); Hg needs q = 0.30 (f'' 6.06, -4.04).
            refine = make_scaling_refine(root, {SILVER: (0.0, 1.218), MERCURY: (0.0, 3.03)})
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine,
                              tolerance=0.3)
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)
            self.assertEqual([(a.selection, a.occupancy_after) for a in result.adjustments],
                             [(MERCURY, 0.3), (SILVER, 0.6)])
            silver_round_2 = next(s for s in result.rounds[1].sites if s.selection == SILVER)
            self.assertEqual(silver_round_2.occupancy, 0.9)   # one site at a time
            final = {outcome.selection: outcome for outcome in result.outcomes}
            self.assertEqual((final[SILVER].occupancy_final, final[MERCURY].occupancy_final), (0.6, 0.3))
            self.assertEqual(occupancy_of(result.final_model, "C", "1", "HG"), 0.3)

    def test_sites_left_out_with_site_are_not_adjusted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)
            refine = make_scaling_refine(root, {SILVER: (0.0, 1.8), MERCURY: (0.0, 3.03)})
            plan = _plan(root, model, mtz, refine=refine, sites=["B:101"])
            result = _execute(plan, root, refine=refine, sites=["B:101"])
            self.assertEqual(result.status, STATUS_CONVERGED)   # Ag within 1 e-
            self.assertEqual(result.adjustments, ())
            self.assertEqual([outcome.selection for outcome in result.outcomes], [SILVER])
            self.assertEqual(occupancy_of(result.final_model, "C", "1", "HG"), 0.5)

    def test_site_whose_fdp_falls_is_held_while_others_are_adjusted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)   # Ag 0.90 and Hg 0.50
            # Hg follows 1/q (theory at q 0.20) only once Ag is below 0.7;
            # before that its f'' falls below q 0.35, as Hg h:1 of PdUHgPdU did.
            refine = make_scaling_refine(root, {
                SILVER: (0.0, 1.218),
                MERCURY: {"S": 2.02, "peak": 0.35, "peak_while": (SILVER, 0.7)},
            })
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine,
                              tolerance=0.3)
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)
            self.assertEqual(
                [(a.round, a.selection, a.occupancy_before, a.occupancy_after, a.method)
                 for a in result.adjustments],
                [(2, MERCURY, 0.5, 0.2, "proportional"),   # f'' 4.04 -> 3.30: fell
                 (3, MERCURY, 0.2, 0.5, "restore"),        # back to the highest f''
                 (4, SILVER, 0.9, 0.6, "proportional"),    # Hg held at 0.50
                 (5, MERCURY, 0.5, 0.2, "proportional")],  # taken up again
            )
            self.assertEqual(
                result.adjustments[1].note,
                "f'' fell from 4.04 to 3.30 as occupancy went 0.50 -> 0.20; back to 0.50 "
                "(highest f'' 4.04), held while the other sites are adjusted",
            )
            self.assertTrue(result.adjustments[3].note.startswith("taken up again after the other sites"))
            mercury_round_4 = next(s for s in result.rounds[3].sites if s.selection == MERCURY)
            self.assertEqual(mercury_round_4.occupancy, 0.5)
            final = {outcome.selection: outcome for outcome in result.outcomes}
            self.assertEqual((final[MERCURY].occupancy_final, final[MERCURY].adjustments,
                              final[MERCURY].held), (0.2, 2, 1))
            self.assertEqual(final[MERCURY].reason,
                             "within +/-0.3 e- of theory after 2 adjustment(s), held back 1 time(s)")
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["rounds"][2]["method"], "restore")
            self.assertEqual(report["occupancy"]["sites"][1]["held"], 1)
            self.assertIn("held there while the other sites are adjusted", report["occupancy"]["method"])

    def test_second_fall_stops_the_site_at_its_best_occupancy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)
            # Hg f'' peaks at 8.66 (q 0.35) whatever happens: theory 10.1 is out of reach.
            refine = make_scaling_refine(root, {
                SILVER: (0.0, 1.218), MERCURY: {"S": 3.03, "peak": 0.35},
            })
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine,
                              tolerance=0.3)
            self.assertEqual(result.status, STATUS_REVIEW)
            self.assertEqual(
                [(a.selection, a.occupancy_after, a.method) for a in result.adjustments],
                [(MERCURY, 0.3, "proportional"),   # f'' 6.06 -> 7.42: rose, fine
                 (MERCURY, 0.17, "secant"),        # 7.42 -> 4.20: fell
                 (MERCURY, 0.3, "restore"),        # held at 0.30
                 (SILVER, 0.6, "proportional"),
                 (MERCURY, 0.22, "proportional"),  # taken up again: 7.42 -> 5.44, fell again
                 (MERCURY, 0.3, "restore")],       # stopped at 0.30
            )
            final = {outcome.selection: outcome for outcome in result.outcomes}
            self.assertEqual(
                final[MERCURY].reason,
                "f'' fell again when occupancy was lowered (0.30 -> 0.22); kept at 0.30, "
                "where f'' was highest (7.42)",
            )
            self.assertEqual((final[MERCURY].occupancy_final, final[MERCURY].held), (0.3, 2))
            self.assertTrue(final[SILVER].converged)
            self.assertEqual(occupancy_of(result.final_model, "C", "1", "HG"), 0.3)
            self.assertEqual(occupancy_of(result.final_model, "B", "101", "AG"), 0.6)
            self.assertEqual(len(result.rounds), 7)

    def test_fdp_rising_as_occupancy_rises_is_held_the_same_way(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)
            refine = make_scaling_refine(root, {
                SILVER: (0.0, 1.8), MERCURY: {"S": 6.06, "trough": 0.52},
            })
            result = _execute(_plan(root, model, mtz, refine=refine), root, refine=refine)
            self.assertEqual(result.status, STATUS_REVIEW)
            self.assertEqual(
                [(a.occupancy_after, a.method) for a in result.adjustments],
                [(0.6, "proportional"), (0.5, "restore"), (0.6, "proportional"), (0.5, "restore")],
            )
            self.assertEqual(
                result.outcomes[1].reason,
                "f'' rose again when occupancy was raised (0.50 -> 0.60); kept at 0.50, "
                "where f'' was lowest (12.12)",
            )

    def test_a_failed_later_round_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_scaling_refine(root, {SILVER: (0.5, 0.6)}, fail_below=0.6)
            result = _execute(_plan(root, model, mtz, refine=refine), root,
                              refine=refine, tolerance=0.1)
            self.assertEqual(result.status, STATUS_FAILED)
            self.assertEqual([item.status for item in result.rounds], [STATUS_ROUND_DONE, STATUS_FAILED])
            self.assertTrue(result.message.startswith("round_002: phenix.refine failed"))
            self.assertIsNone(result.final_model)
            report = json.loads(result.report_path.read_text())
            self.assertEqual((report["status"], report["steps"]["occupancy_adjustment"]),
                             (STATUS_FAILED, "failed"))
            self.assertEqual(len(report["rounds"]), 2)
            self.assertFalse((result.run_directory / "final_model.pdb").exists())


class OccupancyCommandTests(unittest.TestCase):
    _main = AutoFDPCommandTests._main

    def test_command_prints_each_adjustment_and_opens_the_last_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_scaling_refine(root, {SILVER: (0.5, 0.6)})
            code, out, err, launched = self._main(
                root, [str(model), str(mtz), "--tolerance", "0.1"], refine=refine
            )
            self.assertEqual(code, 0, err)
            for text in (
                "AutoFDP step 3: occupancy adjustment (one site at a time, furthest from theory "
                "first; up to 5 adjustment(s) per site",
                "round_002  B:101 AG  occupancy 0.90 -> 0.52 (proportional)  "
                "f'' 1.65 vs 2.03 (-0.38, outside +/-0.1)",
                "round_003  B:101 AG  occupancy 0.52 -> 0.39 (secant)  "
                "f'' 2.04 vs 2.03 (+0.01, within +/-0.1)",
                "Result (round_003)",
                "AUTOFDP_CONVERGED",
            ):
                self.assertIn(text, out)
            row = next(line for line in out.splitlines()
                       if line.startswith("B:101 AG") and "->" in line)
            for value in ("0.90 -> 0.39", "2.04", "+0.01", "after 2 adjustment(s)"):
                self.assertIn(value, row)
            run = mtz.parent.resolve() / "AutoFDP" / "run_001"
            self.assertEqual(launched.call_args.args[0], run / "round_003" / "refined_001.pdb")
            self.assertIn(f"Model with the final occupancies: {run / 'final_model.pdb'}", out)

    def test_command_prints_restore_rounds(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root, mercury=True)
            refine = make_scaling_refine(root, {
                SILVER: (0.0, 1.218), MERCURY: {"S": 3.03, "peak": 0.35},
            })
            code, out, err, _ = self._main(
                root, [str(model), str(mtz), "--tolerance", "0.3", "--no-coot"], refine=refine
            )
            self.assertEqual(code, 2, err)
            self.assertIn("round_004  C:1 HG  occupancy 0.17 -> 0.30 (restore)", out)
            self.assertIn("back to 0.30 (highest f'' 7.42), held while the other sites are adjusted", out)
            self.assertIn("kept at 0.30, where f'' was highest (7.42)", out)

    def test_review_exits_2(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, mtz = make_inputs(root)
            refine = make_scaling_refine(root, {SILVER: (0.0, 2.436)})
            code, out, err, _ = self._main(
                root, [str(model), str(mtz), "--tolerance", "0.1", "--no-coot"], refine=refine
            )
            self.assertEqual(code, 2, err)
            self.assertIn("AUTOFDP_REVIEW", out)
            self.assertIn("f'' above theory at full occupancy", out)


if __name__ == "__main__":
    unittest.main()
