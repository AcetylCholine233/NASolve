import hashlib
import shutil
import tempfile
import unittest
from pathlib import Path

from nasolve.run_context import artifact_reference, dataset_directory, resolve_artifact_path


class RunContextTests(unittest.TestCase):
    def test_autofdp_runs_inside_an_autorefine_round_belong_to_the_dataset(self):
        """AutoFDP from an AutoRefine round writes AutoRefine/round_NNN/AutoFDP/run_NNN
        (user, 2026-11-26); its dataset is the one holding AutoMR/, so the MTZ
        beside it is dataset-anchored, also for an archived PostMR result."""
        with tempfile.TemporaryDirectory() as directory:
            dataset = Path(directory).resolve() / "U-Y_Ag nasolve"
            for nested in (
                dataset / "AutoMR/run_004/AutoRefine/round_006/AutoFDP/run_001",
                dataset / "AutoMR/run_004/PostMRHistory/postmr-002/AutoRefine/round_011/AutoFDP/run_003",
            ):
                nested.mkdir(parents=True)
                self.assertEqual(dataset_directory(nested), dataset)
            mtz = dataset / "staraniso_alldata-unique.mtz"
            mtz.write_bytes(b"MTZ ")
            run = dataset / "AutoMR/run_004/AutoRefine/round_006/AutoFDP/run_001"
            self.assertEqual(artifact_reference(mtz, run),
                             {"anchor": "dataset", "relative_path": "staraniso_alldata-unique.mtz"})
            # Unchanged: AutoFDP beside the MTZ, and metal searches in an AutoMR run.
            self.assertEqual(dataset_directory(dataset / "AutoFDP/run_016"), dataset)
            self.assertEqual(dataset_directory(dataset / "AutoMR/run_004/AutoSolMetal/run_001"), dataset)
            self.assertEqual(dataset_directory(dataset / "AutoMR/run_004"), dataset)

    def test_legacy_collaborator_paths_rebase_to_current_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = Path(directory) / "dataset"
            run = dataset / "AutoMR" / "run_004"
            model = run / "PostMR" / "Model" / "readyset_model.pdb"
            observations = dataset / "staraniso.mtz"
            model.parent.mkdir(parents=True)
            model.write_text("END\n")
            observations.write_bytes(b"mtz")

            old_root = Path("/Users/collaborator/project/dataset")
            self.assertEqual(
                resolve_artifact_path(
                    old_root
                    / "AutoMR"
                    / "run_004"
                    / "PostMR"
                    / "Model"
                    / "readyset_model.pdb",
                    run,
                ),
                model.resolve(),
            )
            self.assertEqual(
                resolve_artifact_path(old_root / "staraniso.mtz", run),
                observations.resolve(),
            )

    def test_rebased_checkout_wins_while_collaborator_path_still_exists(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_run = root / "collaborator" / "dataset" / "AutoMR" / "run_004"
            checkout_run = root / "checkout" / "dataset" / "AutoMR" / "run_004"
            suffix = Path("PostMR") / "Model" / "readyset_model.pdb"
            source_model = source_run / suffix
            checkout_model = checkout_run / suffix
            source_model.parent.mkdir(parents=True)
            checkout_model.parent.mkdir(parents=True)
            source_model.write_text("source\n")
            checkout_model.write_text("checkout\n")

            self.assertEqual(
                resolve_artifact_path(source_model, checkout_run),
                checkout_model.resolve(),
            )

    def test_changed_rebased_file_does_not_fall_back_to_old_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_run = root / "collaborator" / "dataset" / "AutoMR" / "run_004"
            checkout_run = root / "checkout" / "dataset" / "AutoMR" / "run_004"
            suffix = Path("PostMR") / "Model" / "readyset_model.pdb"
            source_model = source_run / suffix
            checkout_model = checkout_run / suffix
            source_model.parent.mkdir(parents=True)
            checkout_model.parent.mkdir(parents=True)
            source_model.write_text("source\n")
            checkout_model.write_text("changed\n")
            source_sha256 = hashlib.sha256(source_model.read_bytes()).hexdigest()

            self.assertIsNone(
                resolve_artifact_path(
                    {"path": str(source_model), "sha256": source_sha256},
                    checkout_run,
                )
            )

    def test_windows_legacy_path_rebases_on_posix_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = Path(directory) / "dataset"
            run = dataset / "AutoMR" / "run_004"
            model = run / "PostMR" / "Model" / "readyset_model.pdb"
            model.parent.mkdir(parents=True)
            model.write_text("END\n")

            self.assertEqual(
                resolve_artifact_path(
                    r"C:\Users\collaborator\NASolve\dataset\AutoMR\run_004"
                    r"\PostMR\Model\readyset_model.pdb",
                    run,
                ),
                model.resolve(),
            )

    def test_new_references_are_anchored_and_cannot_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = Path(directory) / "dataset"
            run = dataset / "AutoMR" / "run_001"
            model = run / "PostMR" / "model.pdb"
            observations = dataset / "data.mtz"
            model.parent.mkdir(parents=True)
            model.write_text("END\n")
            observations.write_bytes(b"mtz")

            self.assertEqual(
                artifact_reference(model, run),
                {"anchor": "run", "relative_path": "PostMR/model.pdb"},
            )
            self.assertEqual(
                artifact_reference(observations, run),
                {"anchor": "dataset", "relative_path": "data.mtz"},
            )
            self.assertEqual(
                resolve_artifact_path(
                    {"anchor": "run", "relative_path": "PostMR/model.pdb"}, run
                ),
                model.resolve(),
            )
            self.assertIsNone(
                resolve_artifact_path(
                    {"anchor": "run", "relative_path": "../../outside.mtz"}, run
                )
            )
            model_ref = artifact_reference(model, run)
            model_ref["sha256"] = hashlib.sha256(model.read_bytes()).hexdigest().upper()
            self.assertEqual(resolve_artifact_path(model_ref, run), model.resolve())

    def test_repository_reference_survives_checkout_relocation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            checkout = root / "checkout"
            run = source / "dataset" / "AutoMR" / "run_004"
            sequence = source / "MR_frames" / "5W6W" / "seq_base.txt"
            run.mkdir(parents=True)
            sequence.parent.mkdir(parents=True)
            (source / "pyproject.toml").write_text("[project]\nname='test'\n")
            sequence.write_text("GAGC\n\nCTGC\n")
            reference = artifact_reference(sequence, run)

            self.assertEqual(
                reference,
                {
                    "anchor": "repository",
                    "relative_path": "MR_frames/5W6W/seq_base.txt",
                },
            )
            shutil.copytree(source, checkout)
            shutil.rmtree(source)
            relocated_run = checkout / "dataset" / "AutoMR" / "run_004"
            self.assertEqual(
                resolve_artifact_path(reference, relocated_run),
                (checkout / "MR_frames" / "5W6W" / "seq_base.txt").resolve(),
            )

    def test_legacy_traversal_is_rejected_even_when_target_exists(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root / "dataset" / "AutoMR" / "run_001"
            run.mkdir(parents=True)
            secret = root / "secret.txt"
            secret.write_text("outside run\n")
            malicious = str(run / ".." / ".." / ".." / "secret.txt")

            self.assertTrue(Path(malicious).is_file())
            self.assertIsNone(resolve_artifact_path(malicious, run))

    def test_ambiguous_legacy_relative_path_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = Path(directory) / "dataset"
            run = dataset / "AutoMR" / "run_001"
            run.mkdir(parents=True)
            (run / "duplicate.mtz").write_bytes(b"run")
            (dataset / "duplicate.mtz").write_bytes(b"dataset")

            self.assertIsNone(resolve_artifact_path("duplicate.mtz", run))


if __name__ == "__main__":
    unittest.main()
