"""``nasolve autorefine RUN --name NAME`` (user's request, 2026-11-26): name a
round's checkpoint as it is made. The folder stays round_NNN (NASolve finds
rounds by number); the name is a bookmark, usable with ``--from NAME``,
``nasolve show RUN --checkpoint NAME`` and ``checkpoints use``, and shown by
``checkpoints list``. Names that look like checkpoint IDs are refused, since a
bookmark is looked up before an ID.

phenix.refine, phenix.mtz.dump and Coot are simulated.
"""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from nasolve.autorefine import AutoRefineError, execute_autorefine
from nasolve.checkpoints import CheckpointError, add_checkpoint
from nasolve.cli import build_parser
from nasolve.coot_view import CootViewError, launch_coot_view

from .test_autorefine import PHENIX_21, make_mtz_dump, make_refine, make_refine_run


class AutoRefineNameTests(unittest.TestCase):
    def refine(self, root: Path, run: Path, **kwargs):
        return execute_autorefine(
            run, make_refine(root), make_mtz_dump(root),
            phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
            macro_cycles=1, auto_select_success=False, **kwargs,
        )

    def test_the_round_gets_the_name_and_later_rounds_can_start_from_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            first = self.refine(root, run, name=" Ag sqp ")
            self.assertEqual(first.status, "AUTOREFINE_READY", first.message)
            self.assertEqual(first.name, "Ag sqp")
            self.assertEqual(first.round_directory.name, "round_001")   # folder unchanged
            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            self.assertEqual(registry["bookmarks"], {"Ag sqp": "refine-001"})
            checkpoint = next(item for item in registry["checkpoints"] if item["id"] == "refine-001")
            self.assertEqual(checkpoint["label"], "Ag sqp")
            self.assertEqual(json.loads(first.report_path.read_text())["name"], "Ag sqp")

            second = self.refine(root, run, from_checkpoint="Ag sqp")
            self.assertEqual(second.parent_checkpoint, "refine-001")
            self.assertIsNone(second.name)
            self.assertIsNone(json.loads(second.report_path.read_text())["name"])

    def test_bad_or_used_names_are_refused_before_anything_is_made(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            self.refine(root, run, name="Ag sqp")
            for name, message in (
                ("Ag sqp", r"Checkpoint name already exists: Ag sqp \(-> refine-001\)"),
                ("refine-007", "looks like a checkpoint ID"),
                ("REFINE-001", "looks like a checkpoint ID"),
                ("postmr", "looks like a checkpoint ID"),
                ("a/b", "may contain letters, numbers, spaces"),
                ("", "must contain 1-80 characters"),
            ):
                with self.subTest(name=name):
                    with self.assertRaisesRegex(AutoRefineError, message):
                        self.refine(root, run, name=name)
            self.assertFalse((run / "AutoRefine" / "round_002").exists())
            # The same rule for checkpoints add --name (it took ID-like names before).
            with self.assertRaisesRegex(CheckpointError, "looks like a checkpoint ID"):
                add_checkpoint(run, name="manual-002")

    def test_show_opens_a_checkpoint_by_its_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            self.refine(root, run, name="Ag sqp")
            captured = {}

            def launcher(command, **kwargs):
                captured["cwd"] = kwargs["cwd"]
                return SimpleNamespace(pid=1)

            result = launch_coot_view(run, Path("/usr/bin/true"), checkpoint="Ag sqp",
                                      launcher=launcher)
            self.assertEqual(result.checkpoint_id, "refine-001")
            self.assertIn("checkpoint refine-001 'Ag sqp'", result.source)
            self.assertEqual(captured["cwd"], run / "CootGUI" / "autorefine" / "refine-001")
            # The ID still works; unknown names and paths do not.
            self.assertEqual(launch_coot_view(run, Path("/usr/bin/true"), checkpoint="refine-001",
                                              launcher=launcher).checkpoint_id, "refine-001")
            with self.assertRaisesRegex(CootViewError, "unknown: no such name"):
                launch_coot_view(run, Path("/usr/bin/true"), checkpoint="no such name",
                                 launcher=launcher)
            with self.assertRaisesRegex(CootViewError, "Not a checkpoint ID or name"):
                launch_coot_view(run, Path("/usr/bin/true"), checkpoint="../refine-001",
                                 launcher=launcher)

    def test_command_line_option(self):
        args = build_parser().parse_args(["autorefine", "run_004", "--name", "Ag sqp"])
        self.assertEqual(args.name, "Ag sqp")
        self.assertIsNone(build_parser().parse_args(["autorefine"]).name)


if __name__ == "__main__":
    unittest.main()
