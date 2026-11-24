"""``nasolve autorefine RUN --model FILE`` (user's request, 2026-11-24).

A round starts from the given PDB instead of the model of the checkpoint it
continues (``--from``, or the current one). Data, phases and restraints still
come from that checkpoint, and it stays the round's parent. The file is
checked like any starting model before anything is created, copied into the
round (``round_NNN/start_model.pdb``, the copy phenix.refine gets) and
recorded with its checksum; the original is not touched.

phenix.refine and phenix.mtz.dump are simulated (the stand-in copies its
input model to refined_001.pdb).
"""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from nasolve.autorefine import AutoRefineError, execute_autorefine

from .helpers import pdb_record
from .test_autorefine import PHENIX_21, make_mtz_dump, make_refine, make_refine_run


def edited_model(path: Path, *, residue: str = "S6G") -> Path:
    path.write_text(
        pdb_record("ATOM", 1, "P", residue, "A", 12, x=1.5)
        + pdb_record("HETATM", 2, "I", "C38", "B", 4, element="I")
        + "END\n"
    )
    return path


class StartModelTests(unittest.TestCase):
    def refine(self, root: Path, run: Path, **kwargs):
        return execute_autorefine(
            run, make_refine(root), make_mtz_dump(root),
            phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
            macro_cycles=1, auto_select_success=False, **kwargs,
        )

    def test_round_starts_from_the_given_model_with_the_checkpoints_data_and_restraints(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            first = self.refine(root, run)
            parent_model = (first.report_path.parent / "refined_001.pdb").resolve()
            edited = edited_model(root / "edited model.pdb")
            original = edited.read_bytes()

            result = self.refine(root, run, from_checkpoint=first.checkpoint_id, start_model=edited)
            round_directory = result.report_path.parent
            copy = (round_directory / "start_model.pdb").resolve()
            self.assertEqual(copy.read_bytes(), original)
            self.assertEqual(edited.read_bytes(), original)          # original untouched
            payload = json.loads(result.report_path.read_text())
            self.assertEqual(payload["parent_checkpoint"], first.checkpoint_id)
            self.assertEqual(payload["inputs"]["model"], str(copy))
            self.assertEqual(payload["inputs"]["start_model"], {
                "source": str(edited.resolve()),
                "sha256": hashlib.sha256(original).hexdigest(),
                "copy": str(copy),
                "instead_of": str(parent_model),
            })
            arguments = json.loads((round_directory / "received_args.json").read_text())
            self.assertIn(str(copy), arguments)
            self.assertNotIn(str(parent_model), arguments)
            self.assertEqual((round_directory / "refined_001.pdb").read_bytes(), original)
            # Same restraints as the checkpoint it continues.
            first_payload = json.loads(first.report_path.read_text())
            self.assertEqual(payload["inputs"]["restraints"], first_payload["inputs"]["restraints"])

            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            checkpoint = next(item for item in registry["checkpoints"] if item["id"] == result.checkpoint_id)
            self.assertEqual(checkpoint["parent"], first.checkpoint_id)
            self.assertEqual(checkpoint["start_model"]["source"], str(edited.resolve()))
            self.assertEqual(checkpoint["start_model"]["relative_path"],
                             f"AutoRefine/{round_directory.name}/start_model.pdb")
            # A round without --model has no start_model record.
            self.assertIsNone(first_payload["inputs"]["start_model"])
            self.assertNotIn("start_model", next(item for item in registry["checkpoints"]
                                                 if item["id"] == first.checkpoint_id))

    def test_unusable_models_are_refused_before_a_round_exists(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            (root / "model.cif").write_text("data_model\n")
            for model, message in (
                (root / "missing.pdb", "Starting model does not exist"),
                (root / "model.cif", r"--model takes a PDB file \(.pdb\), not model.cif"),
                (edited_model(root / "unmutated.pdb", residue="DC"), "A:12 S6G"),
            ):
                with self.subTest(model=model.name), self.assertRaisesRegex(AutoRefineError, message):
                    self.refine(root, run, start_model=model)
            self.assertFalse((run / "AutoRefine" / "round_001").exists())


if __name__ == "__main__":
    unittest.main()
