"""``nasolve autorefine RUN --coordination METAL=drop`` (user's request and
syntax, 2026-11-26): remove a metal's inherited coordination restraints from
this round on, keeping every other restraint. Coordination restraints are
inherited along a branch; before this the only way out was an earlier
checkpoint (--from) with --model.

The dropped metal's restraint file leaves the round's restraint bundle; a metal
that shared that file is written again into this round's metal_coordination.eff
("carried"). With a new METAL=DONOR,... in the same command the metal is
restrained afresh (a replacement).

phenix.refine and phenix.mtz.dump are simulated.
"""

import json
import tempfile
import unittest
from pathlib import Path

from nasolve.autorefine import AutoRefineError, execute_autorefine
from nasolve.cli import build_parser

from .helpers import pdb_record
from .test_autorefine import PHENIX_21, make_mtz_dump, make_refine, make_refine_run


def run_with_two_metals(root: Path) -> Path:
    """Ag Z:1 between C:1 N3 and C:2 O4; Ag Z:2 between C:3 O6 and C:4 N1."""
    run = make_refine_run(root)
    model = run / "PostMR" / "Model" / "readyset_model.pdb"
    kept = model.read_text().replace("END\n", "")
    lines = [
        pdb_record("ATOM", 901, "N3", "DT", "C", 1, element="N", x=-2.2),
        pdb_record("ATOM", 902, "O4", "DT", "C", 2, element="O", x=2.6, y=0.4),
        pdb_record("HETATM", 903, "AG", "AG", "Z", 1, element="AG"),
        pdb_record("ATOM", 904, "O6", "DG", "C", 3, element="O", x=10.0, y=-2.5),
        pdb_record("ATOM", 905, "N1", "DC", "C", 4, element="N", x=10.0, y=2.3),
        pdb_record("HETATM", 906, "AG", "AG", "Z", 2, element="AG", x=10.0),
    ]
    model.write_text(kept + "".join(lines) + "END\n")
    return run


class CoordinationDropTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.run_directory = run_with_two_metals(self.root)

    def refine(self, **kwargs):
        return execute_autorefine(
            self.run_directory, make_refine(self.root), make_mtz_dump(self.root),
            phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
            macro_cycles=1, auto_select_success=False, **kwargs,
        )

    def checkpoint(self, checkpoint_id: str) -> dict:
        registry = json.loads((self.run_directory / "AutoRefine" / "checkpoints.json").read_text())
        return next(item for item in registry["checkpoints"] if item["id"] == checkpoint_id)

    def restraint_paths(self, checkpoint_id: str) -> list[str]:
        return [ref.get("relative_path") for ref in self.checkpoint(checkpoint_id)["restraints"]]

    def test_dropping_the_only_coordination_removes_it_from_here_on(self):
        first = self.refine(coordination=("Z:1=C:1/N3,C:2/O4",))
        first_file = (first.report_path.parent / "metal_coordination.eff").resolve()
        second = self.refine(from_checkpoint=first.checkpoint_id, coordination=("Z:1=drop",))
        self.assertEqual(second.status, "AUTOREFINE_READY", second.message)
        payload = json.loads(second.report_path.read_text())
        self.assertNotIn(str(first_file), payload["command"])
        self.assertFalse((second.report_path.parent / "metal_coordination.eff").exists())
        recorded = payload["inputs"]["metal_coordination"]
        self.assertEqual((recorded["file"], recorded["sites"], recorded["inherited"], recorded["carried"]),
                         (None, [], [], []))
        self.assertEqual([(item["metal"], item["spec"]) for item in recorded["dropped"]],
                         [("AG Z:1 AG", "Z:1=C:1/N3,C:2/O4")])
        self.assertNotIn("metal_coordination", self.checkpoint(second.checkpoint_id))
        self.assertNotIn(f"AutoRefine/{first.report_path.parent.name}/metal_coordination.eff",
                         self.restraint_paths(second.checkpoint_id))
        # Other inherited restraints stay; later rounds inherit no coordination.
        self.assertEqual(len(self.restraint_paths(second.checkpoint_id)),
                         len(self.restraint_paths(first.checkpoint_id)) - 1)
        third = self.refine(from_checkpoint=second.checkpoint_id)
        self.assertEqual(json.loads(third.report_path.read_text())["inputs"]["metal_coordination"]["inherited"], [])
        # The first round's checkpoint is unchanged: a branch from it still has them.
        self.assertEqual([item["metal"] for item in self.checkpoint(first.checkpoint_id)["metal_coordination"]],
                         ["AG Z:1 AG"])

    def test_drop_and_a_new_specification_replace_it(self):
        first = self.refine(coordination=("Z:1=C:1/N3,C:2/O4",))
        second = self.refine(from_checkpoint=first.checkpoint_id,
                             coordination=("Z:1=drop", "Z:1/AG=C:2/O4"))
        self.assertEqual(second.status, "AUTOREFINE_READY", second.message)
        text = (second.report_path.parent / "metal_coordination.eff").read_text()
        self.assertIn("name O4", text)
        self.assertNotIn("name N3", text)
        records = self.checkpoint(second.checkpoint_id)["metal_coordination"]
        self.assertEqual([item["spec"] for item in records], ["Z:1/AG=C:2/O4"])
        self.assertEqual(records[0]["restraint"]["relative_path"],
                         f"AutoRefine/{second.report_path.parent.name}/metal_coordination.eff")

    def test_a_metal_sharing_the_file_is_carried_into_this_round(self):
        first = self.refine(coordination=("Z:1=C:1/N3,C:2/O4", "Z:2=C:3/O6,C:4/N1"))
        first_file = f"AutoRefine/{first.report_path.parent.name}/metal_coordination.eff"
        second = self.refine(from_checkpoint=first.checkpoint_id, coordination=("Z:1=drop",))
        self.assertEqual(second.status, "AUTOREFINE_READY", second.message)
        second_file = second.report_path.parent / "metal_coordination.eff"
        text = second_file.read_text()
        self.assertIn("chain Z and resid 2 and name AG", text)
        self.assertNotIn("chain Z and resid 1 and name AG", text)
        payload = json.loads(second.report_path.read_text())
        self.assertIn(str(second_file.resolve()), payload["command"])
        self.assertNotIn(str((first.report_path.parent / "metal_coordination.eff").resolve()),
                         payload["command"])
        recorded = payload["inputs"]["metal_coordination"]
        self.assertEqual([item["metal"] for item in recorded["carried"]], ["AG Z:2 AG"])
        self.assertEqual(recorded["carried"][0]["carried_over_from"]["relative_path"], first_file)
        self.assertEqual([item["metal"] for item in recorded["dropped"]], ["AG Z:1 AG"])
        records = self.checkpoint(second.checkpoint_id)["metal_coordination"]
        self.assertEqual([(item["metal"], item["restraint"]["relative_path"]) for item in records],
                         [("AG Z:2 AG", f"AutoRefine/{second.report_path.parent.name}/metal_coordination.eff")])
        self.assertNotIn(first_file, self.restraint_paths(second.checkpoint_id))

    def test_what_is_refused_before_a_round_exists(self):
        first = self.refine(coordination=("Z:1=C:1/N3,C:2/O4",))
        rounds = sorted((self.run_directory / "AutoRefine").glob("round_*"))
        for kwargs, message in (
            ({"coordination": ("Z:2=drop",)},
             r"--coordination Z:2=drop: Z:2 has no coordination restraints in checkpoint refine-001 "
             r"\(restrained there: Z:1=C:1/N3,C:2/O4\)"),
            ({"coordination": ("Z:1=drop", "Z:1/AG=DROP")}, "Z:1/AG is dropped twice"),
            ({"coordination": ("Z:1=drop",), "coordination_geometry": ("sqp",)},
             "--geometry needs --coordination"),
            ({"coordination": ("Z:1=C:2/O4",)},
             "already has coordination restraints .* drop them in the same command with "
             "--coordination Z:1=drop"),
            ({"coordination": ("Z1=drop",)}, "is not CHAIN:RESID"),
        ):
            with self.subTest(kwargs=kwargs), self.assertRaisesRegex(AutoRefineError, message):
                self.refine(from_checkpoint=first.checkpoint_id, **kwargs)
        self.assertEqual(sorted((self.run_directory / "AutoRefine").glob("round_*")), rounds)

    def test_command_line_help(self):
        help_text = " ".join(build_parser()._subparsers._group_actions[0].choices["autorefine"]
                             .format_help().split())
        self.assertIn("METAL=drop (e.g. B:101=drop) removes a metal's inherited coordination", help_text)
        args = build_parser().parse_args(["autorefine", "--coordination", "B:101=drop"])
        self.assertEqual(args.coordination, ["B:101=drop"])


if __name__ == "__main__":
    unittest.main()
