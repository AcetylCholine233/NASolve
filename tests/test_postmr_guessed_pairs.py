"""PostMR on a run without a frame: NARestraints guesses every base pair and
restrains them and stacking, with its own defaults (user's request,
2026-11-19: Y-T nasolve AutoMR run_001, nonstandard, got no pair restraints).
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.checkpoints import _initial_restraints
from nasolve.cli import build_parser
from nasolve.postmr import PostMRPreparationError, describe_narestraints, prepare_postmr

from .helpers import make_postmr_report
from .test_autofdp_bstart import bridged_model_text
from .test_postmr import make_data_root, make_ready_set

REPOSITORY = Path(__file__).parents[1]
TRIANGLE = REPOSITORY / "MR_frames" / "5W6W" / "5W6W_noPO4.pdb"   # real 5W6W model, no 5'-phosphate

try:
    HAVE_NARESTRAINTS = importlib.util.find_spec("restraints") is not None
except ImportError:
    HAVE_NARESTRAINTS = False


def frameless_run(root: Path, model: Path, *, frame=None) -> Path:
    run = root / "dataset" / "AutoMR" / "run_001"
    report_path = make_postmr_report(run, model)
    report = json.loads(report_path.read_text())
    report["frame"] = frame
    report["post_mr_plan"] = {"sequences": {}, "standard_pair": None, "mutations": {}}
    report_path.write_text(json.dumps(report))
    return run


def readyset(command, cwd, **kwargs):
    (cwd / "prepared_model.updated.pdb").write_text(Path(command[1]).read_text())
    return SimpleNamespace(returncode=0, stdout="mocked ReadySet")


def prepare(run: Path, root: Path, **kwargs):
    with patch("nasolve.postmr.subprocess.run", side_effect=readyset):
        return prepare_postmr(run, make_ready_set(root), data_root=make_data_root(root), **kwargs)


class GuessedPairWiringTests(unittest.TestCase):
    def test_guessed_restraints_go_to_refinement_and_are_reported(self):
        def builder(compatibility, pair_file, output):
            pair_file.write_text("A 1\nB 1\n")
            output.write_text("refinement.geometry_restraints.edits {}\n")
            return {"retained_pair_count": 1, "retained_pairs": [{"first": "A:1", "second": "B:1"}],
                    "metal_mediated_pairs": ["B:4 DU <-> A:12 DC via AG B:101"]}

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = frameless_run(root, TRIANGLE)
            result = prepare(run, root, guessed_pair_builder=builder)
            phil = result.postmr_directory / "Restraints" / "narestraints_guessed_pairs.phil"
            self.assertIn(phil, result.restraint_paths)
            postmr = json.loads(result.report_path.read_text())
            self.assertEqual(postmr["narestraints"]["mode"], "guessed-pairs")
            # What AutoRefine gives phenix.refine.
            self.assertIn(phil.resolve(), [path.resolve() for path in _initial_restraints(postmr, run)])
            summary = "NARestraints: 1 base pair (guessed from the model) and stacking; " \
                      "1 metal-mediated pair not restrained"
            self.assertEqual(result.narestraints_summary, summary)
            self.assertIn(summary, (result.postmr_directory / "postmr.log").read_text())

    def test_no_narestraints_skips_it_and_is_refused_where_it_does_not_apply(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = frameless_run(root, TRIANGLE)
            called = []
            result = prepare(run, root, no_narestraints=True,
                             guessed_pair_builder=lambda *args: called.append(args))
            self.assertEqual(called, [])
            self.assertEqual(json.loads(result.report_path.read_text())["narestraints"]["mode"], "none")
            self.assertEqual(result.narestraints_summary,
                             "NARestraints: not run; no base-pair or stacking restraints")
        for kwargs, frame in (({"modified_pairs_only": True}, None),
                              ({}, {"name": "W", "aliases": ["W", "5W6W"]})):
            with self.subTest(kwargs=kwargs, frame=frame), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                run = frameless_run(root, TRIANGLE, frame=frame)
                with self.assertRaisesRegex(PostMRPreparationError, "--no-narestraints applies only"):
                    prepare(run, root, no_narestraints=True, **kwargs)
                self.assertFalse((run / "PostMR").exists())

    def test_cli_flag(self):
        self.assertFalse(build_parser().parse_args(["postmr", "run_001"]).no_narestraints)
        self.assertTrue(build_parser().parse_args(["postmr", "run_001", "--no-narestraints"]).no_narestraints)

    def test_summary_for_each_mode(self):
        self.assertEqual(describe_narestraints({"mode": "frame-template", "frame": "W"}),
                         "NARestraints: frame W pair file and secondary-structure restraints, with stacking")
        self.assertEqual(describe_narestraints({"mode": "modified-pairs-only", "retained_pair_count": 2}),
                         "NARestraints: 2 base pairs containing a modified nucleotide (guessed), no stacking")


@unittest.skipUnless(HAVE_NARESTRAINTS, "NARestraints is not installed")
class GuessedPairNARestraintsTests(unittest.TestCase):
    def test_real_triangle_model_gets_pairs_and_stacking(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = frameless_run(root, TRIANGLE)
            result = prepare(run, root)
            record = json.loads(result.report_path.read_text())["narestraints"]
            self.assertEqual(record["mode"], "guessed-pairs")
            pairs = {(pair["first"], pair["second"]) for pair in record["retained_pairs"]}
            self.assertGreaterEqual(len(pairs), 15)   # 19 with NARestraints d1bb52d
            self.assertIn(("A:12", "B:4"), pairs)
            self.assertIn(("A:8", "C:8"), pairs)
            phil = Path(record["restraint_file"]).read_text()
            self.assertIn("stacking_pair", phil)
            self.assertIn("parallelity", phil)

    def test_metal_mediated_pair_is_reported_not_restrained(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "bridged.pdb"
            model.write_text(bridged_model_text())
            from nasolve.postmr import _default_guessed_pair_restraints_builder
            record = _default_guessed_pair_restraints_builder(model, root / "pairs.txt", root / "na.phil")
            self.assertEqual(record["retained_pairs"], [])
            self.assertEqual(len(record["metal_mediated_pairs"]), 1)
            self.assertIn("via AG B:101", record["metal_mediated_pairs"][0])
            self.assertIn("# Metal-mediated pair", (root / "pairs.txt").read_text())


if __name__ == "__main__":
    unittest.main()
