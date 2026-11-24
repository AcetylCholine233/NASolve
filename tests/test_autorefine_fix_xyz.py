"""``nasolve autorefine RUN --fix-xyz CHAIN:RESID[/ATOM]`` (user's request,
2026-11-24): freeze an atom (e.g. a heavy atom) in place for one round.

phenix.refine refines the coordinates of ``refine.sites.individual``; the
fixed atoms are left out of it (``all and not (...)``), which Phenix 2.2.1
honours in coordinate refinement and in global real-space refinement (its
``morph`` step is switched off, and local real-space fitting refits protein
side chains). Occupancy, B and anomalous terms are still refined. The atom's
start and final coordinates are compared; a fixed atom that moved sends the
round to review. Not inherited: give it for each round.

phenix.refine and phenix.mtz.dump are simulated (the stand-in copies its
input model to refined_001.pdb; a wrapper moves the iodine).
"""

import json
import stat
import sys
import tempfile
import unittest
from pathlib import Path

from nasolve.autorefine import AutoRefineError, execute_autorefine

from .test_autorefine import PHENIX_21, make_mtz_dump, make_refine, make_refine_run


def moving_refine(root: Path, shift: float) -> Path:
    """The stand-in phenix.refine, then the iodine moved by ``shift`` A in x."""
    inner = make_refine(root)
    wrapper = root / "phenix.refine.moving"
    wrapper.write_text(
        f"#!{sys.executable}\n"
        "import subprocess, sys\n"
        "from pathlib import Path\n"
        f"code = subprocess.call([{str(inner)!r}, *sys.argv[1:]])\n"
        "model = Path('refined_001.pdb')\n"
        "if model.exists():\n"
        "    lines = model.read_text().splitlines()\n"
        "    out = []\n"
        "    for line in lines:\n"
        "        if line.startswith('HETATM') and line[12:16].strip() == 'I':\n"
        f"            line = line[:30] + '%8.3f' % (float(line[30:38]) + {shift}) + line[38:]\n"
        "        out.append(line)\n"
        "    model.write_text('\\n'.join(out) + '\\n')\n"
        "sys.exit(code)\n"
    )
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)
    return wrapper


class FixXyzTests(unittest.TestCase):
    def refine(self, root: Path, run: Path, refine: Path | None = None, **kwargs):
        return execute_autorefine(
            run, refine or make_refine(root), make_mtz_dump(root),
            phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
            macro_cycles=1, auto_select_success=False, **kwargs,
        )

    def test_fixed_atoms_are_left_out_of_coordinate_refinement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            # B:4 has one atom (I): no atom name needed.
            result = self.refine(root, run, fix_xyz=("B:4", "A:12/P"))
            self.assertEqual(result.status, "AUTOREFINE_READY", result.message)
            params = (result.report_path.parent / "autorefine.params").read_text()
            self.assertIn("      individual = all and not ((chain B and resid 4 and name I) or "
                          "(chain A and resid 12 and name P))\n", params)
            payload = json.loads(result.report_path.read_text())
            self.assertEqual(payload["fixed_xyz"], [
                {"atom": "C38 B:4 I", "selection": "chain B and resid 4 and name I",
                 "start_xyz": [0.0, 0.0, 0.0], "final_xyz": [0.0, 0.0, 0.0], "shift": 0.0},
                {"atom": "S6G A:12 P", "selection": "chain A and resid 12 and name P",
                 "start_xyz": [0.0, 0.0, 0.0], "final_xyz": [0.0, 0.0, 0.0], "shift": 0.0},
            ])
            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            checkpoint = next(item for item in registry["checkpoints"] if item["id"] == result.checkpoint_id)
            self.assertEqual(checkpoint["fixed_xyz"], ["C38 B:4 I", "S6G A:12 P"])

            # Not inherited: the next round refines every atom again.
            later = self.refine(root, run, from_checkpoint=result.checkpoint_id)
            self.assertIn("      individual = all\n",
                          (later.report_path.parent / "autorefine.params").read_text())
            self.assertEqual(json.loads(later.report_path.read_text())["fixed_xyz"], [])

    def test_a_fixed_atom_that_moved_sends_the_round_to_review(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            result = self.refine(root, run, refine=moving_refine(root, 0.3), fix_xyz=("B:4/I",))
            self.assertEqual(result.status, "AUTOREFINE_REVIEW")
            self.assertEqual(result.message, "Refinement completed, but atoms fixed with --fix-xyz "
                                             "moved: C38 B:4 I 0.300 A")
            self.assertEqual(json.loads(result.report_path.read_text())["fixed_xyz"][0]["shift"], 0.3)

    def test_what_is_refused_before_a_round_exists(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            for atoms, message in (
                (("B4",), r"--fix-xyz atom 'B4' is not CHAIN:RESID"),
                (("B:9",), "--fix-xyz atom B:9: no residue B:9 in the model"),
                (("A:12/N1",), r"--fix-xyz atom A:12/N1: no atom N1 in S6G A:12 \(atoms: P\)"),
                (("B:4", "B:4/I"), "--fix-xyz: C38 B:4 I is given twice"),
            ):
                with self.subTest(atoms=atoms), self.assertRaisesRegex(AutoRefineError, message):
                    self.refine(root, run, fix_xyz=atoms)
            self.assertFalse((run / "AutoRefine" / "round_001").exists())


if __name__ == "__main__":
    unittest.main()
