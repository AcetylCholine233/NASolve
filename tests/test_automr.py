import json
import shutil
import tempfile
import unittest
from pathlib import Path

from nasolve.automr import AutoMRInputError, prepare_automr
from nasolve.model_assessment import file_sha256
from nasolve.model_compatibility import load_model_compatibility_facts
from nasolve.run_context import resolve_artifact_path

from .helpers import make_dataset, make_mtz_dump, model_text, w_model_text


VALID = {"1AP", "DT", "DA", "A", "5IU", "DG", "DC", "DF"}
FORCED_W = Path(__file__).parents[1] / "MR_frames/5W6W/5W6W_noPO4.pdb"


class AutoMRPreflightTests(unittest.TestCase):
    def test_mirror_preserves_source_and_freezes_transformed_model(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = make_dataset(Path(directory))

            def fake_mirror(source: Path, destination: Path) -> Path:
                destination.write_text("REMARK MIRRORED\n" + source.read_text())
                return destination

            result = prepare_automr(
                dataset,
                mirror=True,
                mirror_transformer=fake_mirror,
                valid_ligand_codes=VALID,
                mtz_dump_executable=make_mtz_dump(Path(directory) / "tools"),
            )
            source_copy = result.run_directory / "Model" / "source_model.pdb"
            mirrored = result.run_directory / "Model" / "input_model.pdb"
            self.assertEqual(source_copy.read_text(), (dataset / "search.pdb").read_text())
            self.assertTrue(mirrored.read_text().startswith("REMARK MIRRORED"))
            report = json.loads(result.report_path.read_text())
            self.assertTrue(report["inputs"]["mirror"])
            self.assertEqual(report["model_assessment"]["copied_model"], str(mirrored))
            self.assertIn("mirror = true", (result.run_directory / "nasolve.input.txt").read_text())

    def test_mirror_rejects_changed_residue_ids_even_when_counts_match(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = make_dataset(Path(directory))

            def renumbering_mirror(source: Path, destination: Path) -> Path:
                lines = []
                for line in source.read_text().splitlines(keepends=True):
                    if line.startswith(("ATOM  ", "HETATM")) and line[21:22].strip() == "A":
                        resid = int(line[22:26]) + 10
                        line = line[:22] + f"{resid:4d}" + line[26:]
                    lines.append(line)
                destination.write_text("".join(lines))
                return destination

            with self.assertRaisesRegex(AutoMRInputError, "residue IDs by chain"):
                prepare_automr(
                    dataset,
                    mirror=True,
                    mirror_transformer=renumbering_mirror,
                    valid_ligand_codes=VALID,
                    mtz_dump_executable=make_mtz_dump(Path(directory) / "tools"),
                )

    def test_mirror_rejects_lost_atoms(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = make_dataset(Path(directory))

            def lossy_mirror(source: Path, destination: Path) -> Path:
                lines = source.read_text().splitlines(keepends=True)
                destination.write_text("".join(lines[1:]))
                return destination

            with self.assertRaisesRegex(AutoMRInputError, "inventory"):
                prepare_automr(
                    dataset,
                    mirror=True,
                    mirror_transformer=lossy_mirror,
                    valid_ligand_codes=VALID,
                    mtz_dump_executable=make_mtz_dump(Path(directory) / "tools"),
                )

    def test_nonstandard_preflight_generates_input_and_run_records(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = make_dataset(Path(directory))
            result = prepare_automr(
                dataset, valid_ligand_codes=VALID,
                mtz_dump_executable=make_mtz_dump(Path(directory) / "tools"),
            )
            self.assertEqual(result.status, "READY")
            self.assertTrue(result.generated_config)
            self.assertEqual(result.run_directory.name, "run_001")
            self.assertTrue((dataset / "nasolve.txt").is_file())
            for relative in (
                "nasolve.input.txt", "automr.log", "report.json",
                "Model/input_model.pdb", "Model/assessment.json",
                "Model/model_compatibility_facts.json",
            ):
                self.assertTrue((result.run_directory / relative).is_file(), relative)
            report = json.loads(result.report_path.read_text())
            self.assertFalse(report["execution"]["phaser_ran"])
            self.assertTrue(report["model_assessment"]["heteroatoms_preserved"])
            facts = load_model_compatibility_facts(report, result.run_directory)
            self.assertEqual(facts["dimensions"]["frame_identity"]["relation"], "UNKNOWN")
            self.assertEqual(facts["dimensions"]["site_set"]["relation"], "UNKNOWN")
            self.assertEqual(facts["dimensions"]["residue_identity"]["relation"], "UNKNOWN")
            self.assertIsNone(facts["semantics"]["overall_compatibility"])
            self.assertFalse(facts["semantics"]["automatic_reuse_authorized"])

            second = prepare_automr(
                dataset, valid_ligand_codes=VALID,
                mtz_dump_executable=make_mtz_dump(Path(directory) / "tools"),
            )
            self.assertEqual(second.run_directory.name, "run_002")
            self.assertFalse(second.generated_config)

    def test_standard_pair_records_post_mr_plan(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset", include_model=False)
            frames = root / "frames"
            catalogue = frames / "5W6W"
            catalogue.mkdir(parents=True)
            (catalogue / "C_G.pdb").write_text(w_model_text())
            (catalogue / "seq_base.txt").write_text("GAGC\n\nCTGC\n")
            result = prepare_automr(
                dataset, frame_override="W", pair_override="D:T",
                frames_dir=frames, valid_ligand_codes=VALID,
                mtz_dump_executable=make_mtz_dump(root),
            )
            self.assertEqual(result.status, "READY_POST_MR_MUTATION")
            report = json.loads(result.report_path.read_text())
            self.assertEqual(
                report["post_mr_plan"]["standard_pair"]["ligand_codes"],
                ["1AP", "DT"],
            )
            self.assertFalse(
                report["post_mr_plan"]["standard_pair"]["exact_model_match"]
            )
            frozen_sequence = result.run_directory / "Model" / "seq_base.txt"
            self.assertEqual(frozen_sequence.read_text(), "GAGC\n\nCTGC\n")
            sequence_ref = report["inputs"]["frame_sequence"]
            self.assertEqual(
                resolve_artifact_path(sequence_ref, result.run_directory),
                frozen_sequence.resolve(),
            )
            self.assertEqual(sequence_ref["sha256"], file_sha256(frozen_sequence))
            facts = load_model_compatibility_facts(report, result.run_directory)
            self.assertEqual(facts["dimensions"]["frame_identity"]["relation"], "SAME")
            self.assertEqual(facts["dimensions"]["site_set"]["relation"], "UNKNOWN")
            self.assertEqual(
                facts["dimensions"]["symmetry_and_copy_number"]["recipient_symmetry_class"],
                "H3/R3",
            )
            self.assertEqual(
                facts["dimensions"]["terminal_phosphate_chemistry"]["recipient_sites"],
                ["D:1"],
            )

    def test_forced_original_w_scaffold_preserves_recipe_and_sequence_family_delta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset", include_model=False)
            catalogue = root / "frames" / "5W6W"
            catalogue.mkdir(parents=True)
            shutil.copyfile(FORCED_W, catalogue / "5W6W_noPO4.pdb")
            (catalogue / "seq_base.txt").write_text("reviewed-frame-sequence\n")
            (dataset / "nasolve.txt").write_text(
                "[automr]\n"
                "mode = standard\n"
                "frame = W\n"
                "pair = A:T\n"
                "model = 5W6W_noPO4.pdb\n"
                "model_family = w-metal-scaffold\n"
                "sequence_reference = w-metal-scaffold\n"
            )

            result = prepare_automr(
                dataset,
                frames_dir=root / "frames",
                valid_ligand_codes=VALID,
                mtz_dump_executable=make_mtz_dump(root),
            )
            self.assertEqual(result.status, "READY_POST_MR_MUTATION")
            report = json.loads(result.report_path.read_text())

            self.assertEqual(report["inputs"]["model_provider"], {
                "kind": "explicit-standard-model",
                "selection": "user-forced",
                "frame": "W",
                "location": "frame-catalogue",
                "selector": "5W6W_noPO4.pdb",
                "construct_family": "w-metal-scaffold",
            })
            self.assertEqual(report["inputs"]["model_selector"], "5W6W_noPO4.pdb")
            self.assertEqual(report["post_mr_plan"]["allow_op3_sites"], ["D:1"])
            pair = report["post_mr_plan"]["standard_pair"]
            self.assertIsNone(pair["model_pair"])
            self.assertFalse(pair["exact_model_match"])
            self.assertTrue(pair["mutation_required"])
            self.assertEqual(report["sequence_family_preflight"]["differences"], [
                {"site": "A:13", "before": "DC", "after": "DT"},
                {"site": "B:3", "before": "DG", "after": "DA"},
            ])

            input_model = result.run_directory / "Model/input_model.pdb"
            d1 = [
                line for line in input_model.read_text().splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
                and line[21:22].strip() == "D"
                and line[22:26].strip() == "1"
            ]
            self.assertTrue(d1)
            self.assertNotIn("P", {line[12:16].strip() for line in d1})
            self.assertEqual(
                (result.run_directory / "Model/seq_base.txt").read_text(),
                "reviewed-frame-sequence\n",
            )
            snapshot = (result.run_directory / "nasolve.input.txt").read_text()
            self.assertIn("model = 5W6W_noPO4.pdb", snapshot)
            self.assertIn("model_family = w-metal-scaffold", snapshot)
            facts = load_model_compatibility_facts(report, result.run_directory)
            family = facts["dimensions"]["construct_family"]
            self.assertEqual(family["candidate_declared_family"], "w-metal-scaffold")
            self.assertEqual(family["recipient_reference"]["id"], "w-metal-scaffold")
            self.assertEqual(family["relation"], "SAME")
            self.assertIsNone(facts["semantics"]["donor_eligibility"])

    def test_forced_standard_model_missing_required_sites_fails_before_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset", include_model=False)
            catalogue = root / "frames" / "5W6W"
            catalogue.mkdir(parents=True)
            (catalogue / "tiny.pdb").write_text(w_model_text())
            with self.assertRaisesRegex(
                AutoMRInputError, "required PostMR site A:12 is absent"
            ):
                prepare_automr(
                    dataset,
                    frame_override="W",
                    pair_override="D:T",
                    model_override="tiny.pdb",
                    frames_dir=root / "frames",
                    valid_ligand_codes=VALID,
                    mtz_dump_executable=make_mtz_dump(root),
                )
            self.assertFalse((dataset / "AutoMR").exists())

    def test_exact_standard_pair_model_is_ready_without_pair_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset", include_model=False)
            catalogue = root / "frames" / "5W6W"
            catalogue.mkdir(parents=True)
            (catalogue / "C_G.pdb").write_text(w_model_text())
            (catalogue / "A_G.pdb").write_text(w_model_text())
            result = prepare_automr(
                dataset, frame_override="W", pair_override="A:G",
                frames_dir=root / "frames", valid_ligand_codes=VALID,
                mtz_dump_executable=make_mtz_dump(root),
            )
            self.assertEqual(result.status, "READY")
            report = json.loads(result.report_path.read_text())
            self.assertEqual(Path(report["inputs"]["model"]).name, "A_G.pdb")
            self.assertTrue(
                report["post_mr_plan"]["standard_pair"]["exact_model_match"]
            )

    def test_p1_standard_is_blocked_unless_explicitly_shunted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset", include_model=False)
            (dataset / "data_1.cif").write_text(
                "data_test\n_symmetry.space_group_name_H-M 'P 1'\n"
            )
            (dataset / "summary.html").write_text("Spacegroup name P1\n")
            catalogue = root / "frames" / "5W6W"
            catalogue.mkdir(parents=True)
            (catalogue / "C_G.pdb").write_text(w_model_text())
            dump = make_mtz_dump(root, "P 1", 1)
            with self.assertRaisesRegex(AutoMRInputError, "three MR copies"):
                prepare_automr(
                    dataset, frame_override="W", pair_override="D:T",
                    frames_dir=root / "frames", valid_ligand_codes=VALID,
                    mtz_dump_executable=dump,
                )
            result = prepare_automr(
                dataset, frame_override="W", pair_override="D:T",
                frames_dir=root / "frames", valid_ligand_codes=VALID,
                mtz_dump_executable=dump, allow_p1_standard=True,
            )
            self.assertEqual(result.status, "READY_WITH_RED_FLAG")
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["symmetry"]["mr_copies"], 3)
            self.assertIn(
                "allow_p1_standard = true", (dataset / "nasolve.txt").read_text()
            )

    def test_sequence_length_and_mutation_target_are_validated(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = make_dataset(Path(directory))
            config = dataset / "nasolve.txt"
            config.write_text(
                "[automr]\nmode = nonstandard\nmodel = search.pdb\n"
                "[sequences]\nA = A\n"
            )
            with self.assertRaisesRegex(AutoMRInputError, "has length 1"):
                prepare_automr(dataset, valid_ligand_codes=VALID)
            config.write_text(
                "[automr]\nmode = nonstandard\nmodel = search.pdb\n"
                "[mutations]\nA:99 = 5IU\n"
            )
            with self.assertRaisesRegex(AutoMRInputError, "does not exist"):
                prepare_automr(dataset, valid_ligand_codes=VALID)


if __name__ == "__main__":
    unittest.main()
