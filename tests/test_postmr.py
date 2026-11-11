import json
import shutil
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.curated_ligands import validate_curated_dictionary
from nasolve.postmr import (
    MutationAction,
    PostMRPreparationError,
    _coot_script,
    _default_modified_pair_restraints_builder,
    _default_narestraints_builder,
    _filter_scaffold_overlaps,
    _merged_metal_model,
    _modified_nucleotide_sites,
    _patch_narestraints_records,
    _read_report,
    _restore_canonical_mutation_backbones,
    _restore_shared_parent_coordinates,
    _run_readyset,
    _solution_model,
    build_mutation_plan,
    prepare_postmr,
    residue_name,
    scan_anomalous_candidates,
)

from nasolve.model_assessment import file_sha256

from .helpers import make_postmr_report, make_ready_set, pdb_record


PROTECTED_PARENT_ATOMS = (
    "P",
    "OP1",
    "OP2",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C3'",
    "O3'",
    "C2'",
    "C1'",
)


def postmr_model_text(first: str = "DE", second: str = "DG") -> str:
    return "".join([
        pdb_record("HETATM", 1, "P", first, "A", 12),
        pdb_record("ATOM", 2, "P", second, "B", 4),
        pdb_record("HETATM", 3, "MG", "MG", "A", 50, element="MG"),
        "END\n",
    ])


def make_data_root(root: Path) -> Path:
    data = root / "data"
    ligands = data / "ligands"
    restraints = data / "restraints"
    ligands.mkdir(parents=True)
    restraints.mkdir()
    (ligands / "DE.cif").write_text(
        "data_comp_DE\nloop_\n_chem_comp_bond.comp_id\n"
        "_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n"
        "DE C4 S4\nDE P O5'\n#\n"
    )
    package = Path(__file__).parents[1] / "src" / "nasolve" / "data" / "restraints"
    for name in ("5W6W_Std_padd.txt", "5W6W_secondary_structure.eff"):
        shutil.copyfile(package / name, restraints / name)
    return data


class PostMRTests(unittest.TestCase):
    def test_solution_model_prefers_selected_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_run = root / "collaborator" / "dataset" / "AutoMR" / "run_004"
            checkout_run = root / "checkout" / "dataset" / "AutoMR" / "run_004"
            source_model = source_run / "Phaser" / "mr_solution.pdb"
            checkout_model = checkout_run / "Phaser" / "mr_solution.pdb"
            source_model.parent.mkdir(parents=True)
            checkout_model.parent.mkdir(parents=True)
            source_model.write_text("source\n")
            checkout_model.write_text("checkout\n")
            report = {
                "execution": {"phaser": {"solution_pdb": str(source_model)}}
            }

            self.assertEqual(
                _solution_model(checkout_run, report), checkout_model.resolve()
            )

    def test_postmr_starts_from_the_merged_metal_model_when_selected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run_004"
            report_path = make_postmr_report(run, source)
            merged = run / "Phaser" / "mr_solution-merged.pdb"
            merged.write_text(postmr_model_text().replace("END\n", "") + pdb_record(
                "HETATM", 4, "AG", "AG", "Z", 1, element="AG") + "END\n")
            report = json.loads(report_path.read_text())
            report["metal_merge"] = {"schema_version": 1, "selected": merged.name, "models": [{
                "name": merged.name,
                "model": {"anchor": "run", "relative_path": "Phaser/" + merged.name,
                          "sha256": file_sha256(merged)},
            }]}
            report_path.write_text(json.dumps(report))
            seen = []

            def fake_builder(pdb: Path, pairs: Path, output: Path) -> None:
                seen.append(pdb.read_text())
                output.write_text("geometry_restraints.edits {}\n")

            data_root = make_data_root(root)
            # A changed merged model stops before anything is created.
            original_text = merged.read_text()
            merged.write_text(original_text + "REMARK edited\n")
            with self.assertRaisesRegex(PostMRPreparationError, "--original-model"):
                prepare_postmr(run, make_ready_set(root), data_root=data_root,
                               narestraints_builder=fake_builder)
            self.assertFalse((run / "PostMR").exists())
            merged.write_text(original_text)

            result = prepare_postmr(
                run, make_ready_set(root), data_root=data_root,
                narestraints_builder=fake_builder,
            )
            self.assertEqual(result.input_model_source, "metal-merged")
            self.assertEqual(result.input_model, merged.resolve())
            self.assertIn(" AG ", seen[0])
            self.assertIn(" AG ", result.model_path.read_text())
            postmr = json.loads((result.postmr_directory / "report.json").read_text())
            self.assertEqual(postmr["input_model_source"], "metal-merged")
            self.assertEqual(postmr["input_sha256"], file_sha256(merged))
            self.assertEqual(postmr["metal_merge"]["name"], merged.name)

    def test_original_model_option_ignores_a_merged_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run_004"
            report_path = make_postmr_report(run, source)
            report = json.loads(report_path.read_text())
            report["metal_merge"] = {"schema_version": 1, "selected": "gone.pdb", "models": [{
                "name": "gone.pdb", "model": {"anchor": "run", "relative_path": "Phaser/gone.pdb"},
            }]}
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(PostMRPreparationError, "missing or was changed"):
                _merged_metal_model(run, report)
            result = prepare_postmr(
                run, make_ready_set(root), data_root=make_data_root(root),
                narestraints_builder=lambda pdb, pairs, output: output.write_text("x\n"),
                use_merged_model=False,
            )
            self.assertEqual(result.input_model_source, "phaser")
            self.assertEqual(result.input_model.name, "mr_solution.pdb")

    def test_non_object_report_is_rejected_cleanly(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "run_001"
            run.mkdir()
            (run / "report.json").write_text(json.dumps([]))

            with self.assertRaisesRegex(PostMRPreparationError, "not a JSON object"):
                _read_report(run)

    def test_generated_phil_pair_removes_overlapping_project_eff_block(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "project.eff"
            output = root / "filtered.eff"
            source.write_text(
                "refinement {\n  base_pair {\n"
                "    base1 = chain 'A' and resid 12\n"
                "    base2 = chain 'B' and resid 4\n"
                "  }\n  base_pair {\n"
                "    base1 = chain 'A' and resid 3\n"
                "    base2 = chain 'C' and resid 13\n"
                "  }\n}\n"
            )
            report = _filter_scaffold_overlaps(
                source,
                output,
                [{"first": "B:4", "second": "A:12"}],
            )
            filtered = output.read_text()
            self.assertNotIn("resid 12", filtered)
            self.assertIn("resid 13", filtered)
            self.assertEqual(report["removed_overlap_count"], 1)

    def test_default_modified_pair_builder_filters_structured_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = root / "prepared.pdb"
            prepared.write_text("".join([
                pdb_record("HETATM", 1, "C1'", "DE", "A", 12, element="C"),
                pdb_record("HETATM", 2, "C3'", "DE", "A", 12, element="C"),
                pdb_record("ATOM", 3, "C1'", "DC", "B", 4, element="C"),
                pdb_record("ATOM", 4, "C3'", "DC", "B", 4, element="C"),
                "END\n",
            ]))
            compatibility = root / "compatibility.pdb"
            compatibility.write_text(prepared.read_text())
            pair_output = root / "pairs.txt"
            restraint_output = root / "pairs.eff"

            modified = SimpleNamespace(
                first=SimpleNamespace(chain="A", resid="12"),
                second=SimpleNamespace(chain="B", resid="4"),
                recipe="T_C",
                score=0.9,
                noncanonical=True,
            )
            canonical = SimpleNamespace(
                first=SimpleNamespace(chain="C", resid="1"),
                second=SimpleNamespace(chain="D", resid="2"),
                recipe="G_C",
                score=0.8,
                noncanonical=False,
            )
            builder = types.ModuleType("restraints.builder")
            guesser = types.ModuleType("restraints.guesser")
            base_pairs = types.ModuleType("restraints.base_pairs")
            residue_library = types.ModuleType("restraints.residue_library")
            records = [{"Ligand code": "DE", "Base Analog": "T"}]
            builder.load_residue_records = lambda: records
            guesser.load_residue_records = lambda: records
            residue_library.load_residue_records = lambda: records

            def guess_pairs(path: Path, *, allow_noncanonical: bool):
                self.assertEqual(path, compatibility)
                self.assertTrue(allow_noncanonical)
                return [modified, canonical], ["reviewed warning"]

            def write_guess(path: Path, candidates: list[object]):
                self.assertEqual(candidates, [modified])
                path.write_text("A 12\nB 4\n")

            def build_phil_from_pdb(
                path: Path,
                stretches: object,
                output: Path,
                *,
                include_stacking: bool,
            ):
                self.assertEqual(path, compatibility)
                self.assertEqual(stretches, "parsed")
                self.assertFalse(include_stacking)
                output.write_text("geometry_restraints.edits {}\n")

            guesser.guess_pairs = guess_pairs
            guesser.write_guess = write_guess
            builder.build_phil_from_pdb = build_phil_from_pdb
            base_pairs.read_base_pair_file = lambda _path: "parsed"
            package = types.ModuleType("restraints")
            package.builder = builder
            package.guesser = guesser
            modules = {
                "restraints": package,
                "restraints.builder": builder,
                "restraints.guesser": guesser,
                "restraints.base_pairs": base_pairs,
                "restraints.residue_library": residue_library,
            }
            with patch.dict(sys.modules, modules):
                report = _default_modified_pair_restraints_builder(
                    prepared, compatibility, pair_output, restraint_output,
                )
            self.assertEqual(report["guessed_pair_count"], 2)
            self.assertEqual(report["retained_pair_count"], 1)
            self.assertEqual(report["retained_pairs"][0]["first"], "A:12")
            self.assertFalse(report["include_stacking"])
            self.assertTrue(restraint_output.is_file())

    def test_default_narestraints_builder_uses_released_api(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "model.pdb"
            model.write_text(
                pdb_record("ATOM", 1, "P", "DC", "D", 1, element="P") + "END\n"
            )
            pairs = root / "pairs.txt"
            pairs.write_text("D 1\nD 1\n")
            output = root / "restraints.phil"
            seen: dict[str, object] = {}

            builder = types.ModuleType("restraints.builder")
            base_pairs = types.ModuleType("restraints.base_pairs")
            residue_library = types.ModuleType("restraints.residue_library")
            builder.load_residue_records = lambda: []
            residue_library.load_residue_records = lambda: []
            base_pairs.read_base_pair_file = lambda _path: "parsed"

            def build_phil_from_pdb(
                path: Path, stretches: object, destination: Path, *,
                include_stacking: bool,
            ) -> None:
                seen["path"] = path
                seen["stretches"] = stretches
                seen["include_stacking"] = include_stacking
                destination.write_text("geometry_restraints.edits {}\n")

            builder.build_phil_from_pdb = build_phil_from_pdb
            package = types.ModuleType("restraints")
            package.builder = builder
            modules = {
                "restraints": package,
                "restraints.builder": builder,
                "restraints.base_pairs": base_pairs,
                "restraints.residue_library": residue_library,
            }
            with patch.dict(sys.modules, modules):
                corrections = _default_narestraints_builder(
                    model, pairs, output
                )

            self.assertEqual(corrections, [])
            self.assertEqual(seen["path"], model)
            self.assertEqual(seen["stretches"], "parsed")
            self.assertTrue(seen["include_stacking"])
            self.assertTrue(output.is_file())

    def test_mirrored_canonical_targets_do_not_revert_to_d_dna(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "model.pdb"
            model.write_text(postmr_model_text("0DC", "0DG"))
            run = root / "run"
            report_path = make_postmr_report(
                run, model, first="DC", second="DG", requested="C:G",
            )
            report = json.loads(report_path.read_text())
            report["inputs"] = {"mirror": True}
            actions = build_mutation_plan(report, model)
            self.assertEqual(
                [(action.after, action.method) for action in actions],
                [("0DC", "none"), ("0DG", "none")],
            )

            report["post_mr_plan"]["standard_pair"]["ligand_codes"][0] = "DA"
            with self.assertRaisesRegex(
                PostMRPreparationError, "unmirror/Coot/remirror"
            ):
                build_mutation_plan(report, model)

    def test_modified_site_detection_uses_chemical_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / "model.pdb"
            model.write_text("".join([
                pdb_record("HETATM", 1, "C1'", "DE", "A", 12, element="C"),
                pdb_record("HETATM", 2, "C3'", "DE", "A", 12, element="C"),
                pdb_record("ATOM", 3, "C1'", "DC", "B", 4, element="C"),
                pdb_record("ATOM", 4, "C3'", "DC", "B", 4, element="C"),
                "END\n",
            ]))
            self.assertEqual(_modified_nucleotide_sites(model), {"A:12"})

    def test_readyset_successful_noop_preserves_and_audits_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "prepared_model.pdb"
            model.write_text(postmr_model_text("DC", "DG"))
            readyset_directory = root / "ReadySet"
            readyset_directory.mkdir()
            executable = root / "phenix.ready_set"
            executable.write_text(
                "#!/bin/sh\n"
                "echo 'No unknown residues'\n"
                "exit 0\n"
            )
            executable.chmod(executable.stat().st_mode | 0o100)

            checked, log, generated_cif, command, audit = _run_readyset(
                model, executable, readyset_directory, None, None
            )

            self.assertTrue(checked.is_file())
            self.assertEqual(checked.read_text(), model.read_text())
            self.assertIn("No unknown residues", log.read_text())
            self.assertIsNone(generated_cif)
            self.assertEqual(command[0], str(executable))
            self.assertEqual(audit["readyset_output_mode"], "successful-noop")
            self.assertFalse((readyset_directory / "prepared_model.updated.pdb").exists())

    def test_modified_pairs_only_works_for_nonstandard_and_still_runs_readyset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text("".join([
                pdb_record("HETATM", 1, "P", "DE", "A", 12, element="P"),
                pdb_record("HETATM", 2, "C1'", "DE", "A", 12, element="C"),
                pdb_record("HETATM", 3, "C3'", "DE", "A", 12, element="C"),
                pdb_record("ATOM", 4, "P", "DC", "B", 4, element="P"),
                pdb_record("ATOM", 5, "C1'", "DC", "B", 4, element="C"),
                pdb_record("ATOM", 6, "C3'", "DC", "B", 4, element="C"),
                "END\n",
            ]))
            run = root / "run"
            report_path = make_postmr_report(run, source)
            report = json.loads(report_path.read_text())
            report["frame"] = None
            report["post_mr_plan"] = {
                "sequences": {}, "standard_pair": None, "mutations": {},
            }
            report_path.write_text(json.dumps(report))

            def fake_modified_builder(
                prepared: Path,
                compatibility: Path,
                pairs: Path,
                output: Path,
            ) -> dict[str, object]:
                self.assertIn(" DE ", prepared.read_text())
                self.assertIn(" DE ", compatibility.read_text())
                pairs.write_text("A 12\nB 4\n")
                output.write_text("geometry_restraints.edits {}\n")
                return {
                    "mode": "modified-pairs-only",
                    "modified_sites": ["A:12"],
                    "retained_pair_count": 1,
                    "retained_pairs": [{"first": "A:12", "second": "B:4"}],
                    "compatibility_corrections": [],
                }

            result = prepare_postmr(
                run,
                make_ready_set(root),
                data_root=make_data_root(root),
                modified_pairs_only=True,
                modified_pair_builder=fake_modified_builder,
            )
            restraint_names = {path.name for path in result.restraint_paths}
            self.assertIn("narestraints_modified_pairs.phil", restraint_names)
            self.assertIn("DE.cif", restraint_names)
            self.assertNotIn("5W6W_secondary_structure.eff", restraint_names)
            self.assertTrue(result.readyset_log.is_file())
            postmr_report = json.loads(result.report_path.read_text())
            self.assertEqual(
                postmr_report["narestraints"]["mode"], "modified-pairs-only"
            )
            self.assertEqual(
                postmr_report["narestraints"]["retained_pair_count"], 1
            )

    def test_modified_pairs_only_empty_selection_is_successful_noop(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run"
            report_path = make_postmr_report(run, source)
            report = json.loads(report_path.read_text())
            report["frame"] = None
            report["post_mr_plan"] = {
                "sequences": {}, "standard_pair": None, "mutations": {},
            }
            report_path.write_text(json.dumps(report))

            def no_pairs(*args: Path) -> dict[str, object]:
                args[2].write_text("")
                return {
                    "mode": "modified-pairs-only",
                    "modified_sites": [],
                    "retained_pair_count": 0,
                    "retained_pairs": [],
                    "compatibility_corrections": [],
                }

            result = prepare_postmr(
                run,
                make_ready_set(root),
                data_root=make_data_root(root),
                modified_pairs_only=True,
                modified_pair_builder=no_pairs,
            )
            self.assertFalse(
                (result.postmr_directory / "Restraints" / "narestraints_modified_pairs.phil").exists()
            )
            self.assertTrue(result.readyset_log.is_file())

    def test_modified_pair_phil_and_project_eff_compose_for_w_frame(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run"
            make_postmr_report(run, source)

            def one_pair(
                prepared: Path,
                compatibility: Path,
                pairs: Path,
                output: Path,
            ) -> dict[str, object]:
                pairs.write_text("A 12\nB 4\n")
                output.write_text("geometry_restraints.edits {}\n")
                return {
                    "mode": "modified-pairs-only",
                    "modified_sites": ["A:12"],
                    "retained_pair_count": 1,
                    "retained_pairs": [{"first": "A:12", "second": "B:4"}],
                    "compatibility_corrections": [],
                }

            result = prepare_postmr(
                run,
                make_ready_set(root),
                data_root=make_data_root(root),
                modified_pairs_only=True,
                modified_pair_builder=one_pair,
            )
            names = {path.name for path in result.restraint_paths}
            self.assertIn("narestraints_modified_pairs.phil", names)
            self.assertIn("5W6W_secondary_structure.eff", names)
            report = json.loads(result.report_path.read_text())
            self.assertTrue(
                report["narestraints"]["project_scaffold"]
                ["generated_phil_is_authoritative"]
            )

    def test_full_dna_sequence_expands_by_frozen_chain_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "model.pdb"
            model.write_text("".join([
                pdb_record("ATOM", 1, "P", "DA", "A", 5),
                pdb_record("ATOM", 2, "P", "DC", "A", 7),
                "END\n",
            ]))
            run = root / "run"
            report_path = make_postmr_report(run, model, first="DA", second="DC")
            report = json.loads(report_path.read_text())
            report["post_mr_plan"] = {
                "sequences": {"A": "GT"},
                "standard_pair": None,
                "mutations": {"A:7": {"requested": "iU", "ligand_code": "5IU"}},
            }
            report["model_assessment"] = {
                "polymer_residue_ids_by_chain": {"A": ["5", "7"]},
            }
            actions = build_mutation_plan(report, model)
            self.assertEqual(
                [(action.site, action.before, action.after, action.method) for action in actions],
                [
                    ("A:5", "DA", "DG", "coot-mutate-base"),
                    ("A:7", "DC", "5IU", "coot-parent-overlap"),
                ],
            )
            self.assertEqual(actions[1].parent_code, "DT")
            script = _coot_script(
                model,
                root / "after.pdb",
                actions,
                {"5IU": root / "5IU.cif"},
                {"A:7": root / "parent_A_7_DT.pdb"},
            )
            self.assertIn("mutate_base(imol, 'A', 5, '', 'DG')", script)
            self.assertIn("mutate_base(imol, 'A', 7, '', 'DT')", script)

    def test_canonical_sequence_mutation_restores_original_backbone(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            mutated = root / "mutated.pdb"
            output = root / "output.pdb"
            source.write_text("".join([
                pdb_record(
                    "ATOM", 1, "P", "DA", "A", 5, element="P",
                    x=1.0, y=2.0, z=3.0, occupancy=0.75, b_factor=17.0,
                ),
                pdb_record(
                    "ATOM", 2, "O1P", "DA", "A", 5, element="O",
                    x=2.0, y=3.0, z=4.0, occupancy=0.75, b_factor=17.0,
                ),
                pdb_record(
                    "ATOM", 3, "C1'", "DA", "A", 5, element="C",
                    x=3.0, y=4.0, z=5.0, occupancy=0.75, b_factor=17.0,
                ),
                "END\n",
            ]))
            mutated.write_text("".join([
                pdb_record("ATOM", 1, "P", "DG", "A", 5, element="P", x=101.0),
                pdb_record("ATOM", 2, "OP1", "DG", "A", 5, element="O", x=102.0),
                pdb_record("ATOM", 3, "C1'", "DG", "A", 5, element="C", x=103.0),
                pdb_record("ATOM", 4, "N9", "DG", "A", 5, element="N", x=104.0),
                "END\n",
            ]))
            action = MutationAction("A:5", "DA", "DG", "coot-mutate-base")
            report = _restore_canonical_mutation_backbones(
                source, mutated, output, (action,)
            )
            records = {
                line[12:16].strip(): line
                for line in output.read_text().splitlines()
                if line.startswith(("ATOM", "HETATM"))
            }
            self.assertEqual(float(records["P"][30:38]), 1.0)
            self.assertEqual(float(records["OP1"][30:38]), 2.0)
            self.assertEqual(float(records["C1'"][30:38]), 3.0)
            self.assertEqual(float(records["N9"][30:38]), 104.0)
            self.assertEqual(report["A:5"]["source"], "pre-Coot input model")

    def test_w_frame_keeps_de_and_runs_readyset_without_coot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run_004"
            make_postmr_report(run, source)
            calls = []

            def fake_builder(pdb: Path, pairs: Path, output: Path) -> None:
                calls.append((pdb.read_text(), pairs.read_text()))
                output.write_text("geometry_restraints.edits {}\n")

            result = prepare_postmr(
                run,
                make_ready_set(root),
                data_root=make_data_root(root),
                narestraints_builder=fake_builder,
            )
            self.assertEqual(result.status, "POSTMR_READY")
            self.assertEqual(residue_name(result.model_path, "A:12"), "DE")
            self.assertEqual(residue_name(result.model_path, "B:4"), "DG")
            self.assertIn(" DE ", calls[0][0])
            self.assertEqual(calls[0][1], "A 11:13\nB 5:3\n")
            self.assertFalse((result.postmr_directory / "Coot" / "coot.log").exists())
            readyset = (result.postmr_directory / "ReadySet" / "ready_set.log").read_text()
            self.assertIn("user provided restraints", readyset)
            report = json.loads((run / "report.json").read_text())
            self.assertEqual(report["stage"], "postmr")
            self.assertEqual(report["postmr"]["mutation_actions"][0]["method"], "none")
            self.assertIsNone(report["postmr"]["component_identity"]["DE"]["deposition_code"])

    def test_secondary_template_has_base_pairs_but_no_stacking_pairs(self):
        template = (
            Path(__file__).parents[1]
            / "src" / "nasolve" / "data" / "restraints"
            / "5W6W_secondary_structure.eff"
        ).read_text()
        self.assertEqual(template.count("base_pair {"), 17)
        self.assertNotIn("stacking_pair", template)
        self.assertNotIn("file_name", template)
        self.assertNotIn("unit_cell", template)

    def test_legacy_e_8ro_report_is_migrated_to_de(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run_legacy"
            make_postmr_report(run, source, first="8RO")

            def fake_builder(_pdb: Path, _pairs: Path, output: Path) -> None:
                output.write_text("geometry_restraints.edits {}\n")

            result = prepare_postmr(
                run,
                make_ready_set(root),
                data_root=make_data_root(root),
                narestraints_builder=fake_builder,
            )
            self.assertEqual(residue_name(result.model_path, "A:12"), "DE")

    def test_fd_pair_uses_canonical_parents_and_dictionary_overlap(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "model.pdb"
            model.write_text(postmr_model_text("DC", "DG"))
            run = root / "run"
            report_path = make_postmr_report(
                run,
                model,
                first="DF",
                second="1AP",
                requested="F:D",
            )
            report = json.loads(report_path.read_text())
            actions = build_mutation_plan(report, model)
            self.assertEqual(
                [
                    (action.site, action.after, action.method, action.parent_code)
                    for action in actions
                ],
                [
                    ("A:12", "DF", "coot-parent-overlap", "DT"),
                    ("B:4", "1AP", "coot-parent-overlap", "DA"),
                ],
            )
            self.assertEqual(actions[0].deposition_code, "A1AAZ")
            self.assertEqual(actions[1].deposition_code, "1AP")

            script = _coot_script(
                model,
                root / "raw.pdb",
                actions,
                {"DF": root / "DF.cif", "1AP": root / "1AP.cif"},
                {"A:12": root / "A12-DT.pdb", "B:4": root / "B4-DA.pdb"},
            )
            self.assertIn("mutate_base(imol, 'A', 12, '', 'DT')", script)
            self.assertIn("mutate_base(imol, 'B', 4, '', 'DA')", script)
            self.assertIn("get_monomer_from_dictionary('DF', 0)", script)
            self.assertIn("get_monomer_from_dictionary('1AP', 0)", script)
            self.assertIn("overlap_ligands_py", script)
            self.assertIn("add_ligand_delete_residue_copy_molecule", script)
            self.assertIn("ligand_imol, 'A', 1", script)
            self.assertIn("delete_hydrogen_atoms", script)

    def test_shared_parent_coordinates_are_restored(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            parent = root / "parent.pdb"
            raw = root / "raw.pdb"
            output = root / "restored.pdb"
            parent_lines = []
            raw_lines = []
            for serial, atom in enumerate(PROTECTED_PARENT_ATOMS, 1):
                element = "P" if atom == "P" else atom[0]
                raw_atom = {"OP1": "O1P", "OP2": "O2P"}.get(atom, atom)
                parent_lines.append(
                    pdb_record(
                        "ATOM",
                        serial,
                        atom,
                        "DT",
                        "A",
                        12,
                        element=element,
                        x=float(serial),
                        y=float(serial + 1),
                        z=float(serial + 2),
                        occupancy=0.75,
                        b_factor=17.0,
                    )
                )
                raw_lines.append(
                    pdb_record(
                        "HETATM",
                        serial,
                        raw_atom,
                        "DF",
                        "A",
                        12,
                        element=element,
                        x=float(serial + 100),
                        y=float(serial + 101),
                        z=float(serial + 102),
                    )
                )
            parent_lines.extend([
                pdb_record(
                    "ATOM", 90, "C2", "DT", "A", 12, element="C", x=10.0,
                    occupancy=0.75, b_factor=17.0,
                ),
                pdb_record(
                    "ATOM", 91, "O2", "DT", "A", 12, element="O", x=11.0,
                    occupancy=0.50, b_factor=33.0,
                ),
            ])
            raw_lines.extend([
                pdb_record(
                    "HETATM", 90, "C2", "DF", "A", 12, element="C",
                    x=100.0, y=100.0,
                ),
                pdb_record(
                    "HETATM", 99, "S1", "DF", "A", 12, element="S",
                    x=100.0, y=102.0,
                ),
            ])
            parent.write_text("".join(parent_lines) + "END\n")
            raw.write_text("".join(raw_lines) + "END\n")
            action = MutationAction(
                "A:12",
                "DC",
                "DF",
                "coot-parent-overlap",
                parent_code="DT",
                deposition_code="A1AAZ",
            )
            restoration = _restore_shared_parent_coordinates(
                raw,
                output,
                (action,),
                {"A:12": parent},
            )
            records = {
                line[12:16].strip(): line
                for line in output.read_text().splitlines()
                if line.startswith(("ATOM", "HETATM"))
            }
            self.assertEqual(float(records["P"][30:38]), 1.0)
            self.assertEqual(float(records["P"][54:60]), 0.75)
            self.assertEqual(float(records["P"][60:66]), 17.0)
            self.assertEqual(float(records["O1P"][30:38]), 2.0)
            self.assertEqual(float(records["O2P"][30:38]), 3.0)
            self.assertEqual(float(records["S1"][30:38]), 12.0)
            self.assertEqual(float(records["S1"][38:46]), 0.0)
            self.assertEqual(float(records["S1"][54:60]), 0.50)
            self.assertEqual(float(records["S1"][60:66]), 33.0)
            self.assertEqual(restoration["A:12"]["count"], 12)
            self.assertEqual(
                restoration["A:12"]["substitutions"],
                [{
                    "parent_atom": "O2",
                    "target_atom": "S1",
                    "anchor_atom": "C2",
                    "bond_length": 2.0,
                    "placement": "canonical-parent-vector",
                }],
            )

    def test_curated_topology_validation_rejects_missing_or_forbidden_bonds(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            missing = root / "missing.cif"
            missing.write_text(
                "data_DF\nloop_\n_chem_comp_bond.atom_id_1\n"
                "_chem_comp_bond.atom_id_2\nC2 N3\n#\n"
            )
            with self.assertRaisesRegex(ValueError, "missing C2-S1"):
                validate_curated_dictionary("DF", missing)

            forbidden = root / "forbidden.cif"
            forbidden.write_text(
                "data_DF\nloop_\n_chem_comp_bond.atom_id_1\n"
                "_chem_comp_bond.atom_id_2\nC2 S1\nN3 S1\n#\n"
            )
            with self.assertRaisesRegex(ValueError, "forbidden N3-S1"):
                validate_curated_dictionary("DF", forbidden)

    def test_s6g_sulfur_is_projected_along_parent_o6_vector(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            parent = root / "parent-dg.pdb"
            raw = root / "raw-s6g.pdb"
            output = root / "restored-s6g.pdb"
            parent_lines = []
            raw_lines = []
            for serial, atom in enumerate(PROTECTED_PARENT_ATOMS, 1):
                element = "P" if atom == "P" else atom[0]
                parent_lines.append(
                    pdb_record("ATOM", serial, atom, "DG", "B", 4, element=element)
                )
                raw_lines.append(
                    pdb_record(
                        "HETATM", serial, atom, "S6G", "B", 4,
                        element=element, x=100.0,
                    )
                )
            parent_lines.extend([
                pdb_record("ATOM", 90, "C6", "DG", "B", 4, element="C", x=5.0),
                pdb_record(
                    "ATOM", 91, "O6", "DG", "B", 4, element="O", x=5.0, y=1.0,
                    occupancy=0.60, b_factor=24.0,
                ),
            ])
            raw_lines.extend([
                pdb_record(
                    "HETATM", 90, "C6", "S6G", "B", 4,
                    element="C", x=100.0, y=100.0,
                ),
                pdb_record(
                    "HETATM", 91, "S6", "S6G", "B", 4,
                    element="S", x=98.0, y=100.0,
                ),
            ])
            parent.write_text("".join(parent_lines) + "END\n")
            raw.write_text("".join(raw_lines) + "END\n")
            action = MutationAction(
                "B:4", "DG", "S6G", "coot-parent-overlap",
                parent_code="DG", deposition_code="S6G",
            )
            restoration = _restore_shared_parent_coordinates(
                raw, output, (action,), {"B:4": parent}
            )
            records = {
                line[12:16].strip(): line
                for line in output.read_text().splitlines()
                if line.startswith(("ATOM", "HETATM"))
            }
            self.assertEqual(float(records["S6"][30:38]), 5.0)
            self.assertEqual(float(records["S6"][38:46]), 2.0)
            self.assertEqual(float(records["S6"][54:60]), 0.60)
            self.assertEqual(
                restoration["B:4"]["substitutions"][0]["placement"],
                "canonical-parent-vector",
            )

    def test_legacy_df_narestraints_atom_mapping_is_patched_process_locally(self):
        original = [
            {"Ligand code": "DF", "Base Analog": "T", "O2": "S2"},
            {"Ligand code": "1AP", "Base Analog": "D", "N2": "N2"},
        ]
        patched, corrections = _patch_narestraints_records(original, {"DF", "1AP"})
        self.assertEqual(original[0]["O2"], "S2")
        self.assertEqual(patched[0]["O2"], "S1")
        self.assertEqual(
            corrections,
            [{
                "ligand_code": "DF",
                "canonical_atom": "O2",
                "before": "S2",
                "after": "S1",
                "scope": "process-local",
            }],
        )

        already_fixed, corrections = _patch_narestraints_records(
            [{"Ligand code": "DF", "O2": "S1"}], {"DF"}
        )
        self.assertEqual(already_fixed[0]["O2"], "S1")
        self.assertEqual(corrections, [])

        s6g, corrections = _patch_narestraints_records(
            [{"Ligand code": "S6G", "C5": "C5 "}], {"S6G"}
        )
        self.assertEqual(s6g[0]["C5"], "C5")
        self.assertEqual(corrections[0]["before"], "C5 ")
        self.assertEqual(corrections[0]["after"], "C5")

    def test_q_uses_s6g_dictionary_overlap_from_dg_parent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "model.pdb"
            model.write_text(postmr_model_text("DC", "DG"))
            run = root / "run"
            report_path = make_postmr_report(
                run,
                model,
                first="DC",
                second="S6G",
                requested="C:Q",
            )
            report = json.loads(report_path.read_text())
            actions = build_mutation_plan(report, model)
            self.assertEqual(actions[0].method, "none")
            self.assertEqual(
                (actions[1].after, actions[1].method, actions[1].parent_code),
                ("S6G", "coot-parent-overlap", "DG"),
            )

    def test_q_ic_pair_uses_dg_and_dc_parents(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "model.pdb"
            model.write_text(postmr_model_text("DC", "DG"))
            run = root / "run"
            report_path = make_postmr_report(
                run,
                model,
                first="S6G",
                second="C38",
                requested="Q:iC",
            )
            actions = build_mutation_plan(json.loads(report_path.read_text()), model)
            self.assertEqual(
                [(action.after, action.parent_code) for action in actions],
                [("S6G", "DG"), ("C38", "DC")],
            )

    def test_c38_iodine_is_projected_outward_from_c5(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            parent = root / "parent-dc.pdb"
            raw = root / "raw-c38.pdb"
            output = root / "restored-c38.pdb"
            dictionary = root / "C38.cif"
            dictionary.write_text(
                "data_C38\nloop_\n"
                "_chem_comp_atom.comp_id\n"
                "_chem_comp_atom.atom_id\n"
                "_chem_comp_atom.pdbx_model_Cartn_x_ideal\n"
                "_chem_comp_atom.pdbx_model_Cartn_y_ideal\n"
                "_chem_comp_atom.pdbx_model_Cartn_z_ideal\n"
                "C38 C5 0.0 0.0 0.0\n"
                "C38 I 0.0 2.095 0.0\n#\n"
            )
            parent_lines = []
            raw_lines = []
            for serial, atom in enumerate(PROTECTED_PARENT_ATOMS, 1):
                element = "P" if atom == "P" else atom[0]
                parent_lines.append(
                    pdb_record("ATOM", serial, atom, "DC", "B", 4, element=element)
                )
                raw_lines.append(
                    pdb_record(
                        "HETATM", serial, atom, "C38", "B", 4,
                        element=element, x=100.0,
                    )
                )
            parent_lines.extend([
                pdb_record("ATOM", 90, "C4", "DC", "B", 4, element="C", x=-1.0),
                pdb_record(
                    "ATOM", 91, "C5", "DC", "B", 4, element="C", y=1.0,
                    occupancy=0.65, b_factor=27.0,
                ),
                pdb_record("ATOM", 92, "C6", "DC", "B", 4, element="C", x=1.0),
            ])
            raw_lines.extend([
                pdb_record("HETATM", 90, "C4", "C38", "B", 4, element="C"),
                pdb_record(
                    "HETATM", 91, "C5", "C38", "B", 4,
                    element="C", x=100.0, y=100.0,
                ),
                pdb_record(
                    "HETATM", 93, "I", "C38", "B", 4,
                    element="I", x=101.514, y=100.0,
                ),
                pdb_record("HETATM", 92, "C6", "C38", "B", 4, element="C"),
            ])
            parent.write_text("".join(parent_lines) + "END\n")
            raw.write_text("".join(raw_lines) + "END\n")
            action = MutationAction(
                "B:4", "DG", "C38", "coot-parent-overlap",
                parent_code="DC", deposition_code="C38",
            )
            restoration = _restore_shared_parent_coordinates(
                raw,
                output,
                (action,),
                {"B:4": parent},
                {"C38": dictionary},
            )
            records = {
                line[12:16].strip(): line
                for line in output.read_text().splitlines()
                if line.startswith(("ATOM", "HETATM"))
            }
            self.assertAlmostEqual(float(records["I"][30:38]), 0.0, places=3)
            self.assertAlmostEqual(float(records["I"][38:46]), 3.095, places=3)
            self.assertEqual(float(records["I"][54:60]), 0.65)
            self.assertEqual(float(records["I"][60:66]), 27.0)
            self.assertEqual(
                restoration["B:4"]["substitutions"][0]["placement"],
                "canonical-ring-outward-bisector",
            )
            self.assertEqual(
                restoration["B:4"]["substitutions"][0]["bond_length_source"],
                "dictionary-ideal-coordinates",
            )

    def test_anomalous_scan_is_element_driven_and_nucleotide_scoped(self):
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / "heavy-atoms.pdb"
            model.write_text("".join([
                pdb_record("HETATM", 1, "C1'", "C38", "A", 12, element="C"),
                pdb_record("HETATM", 2, "O4'", "C38", "A", 12, element="O"),
                pdb_record("HETATM", 3, "I", "C38", "A", 12, element="I"),
                pdb_record("HETATM", 4, "C1'", "NEW", "B", 4, element="C"),
                pdb_record("HETATM", 5, "C3'", "NEW", "B", 4, element="C"),
                pdb_record("HETATM", 6, "BR", "NEW", "B", 4, element=""),
                pdb_record("HETATM", 7, "I", "IOD", "C", 1, element="I"),
                pdb_record("HETATM", 8, "S6", "S6G", "D", 4, element="S"),
                "END\n",
            ]))
            candidates = scan_anomalous_candidates(model)
            self.assertEqual(
                [(item["site"], item["element"], item["element_source"]) for item in candidates],
                [
                    ("A:12", "I", "pdb-element-column"),
                    ("B:4", "BR", "atom-name-fallback"),
                ],
            )

    def test_review_result_requires_explicit_permission(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run"
            report_path = make_postmr_report(run, source)
            report = json.loads(report_path.read_text())
            report["status"] = "MR_REVIEW"
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(PostMRPreparationError, "--allow-mr-review"):
                prepare_postmr(
                    run,
                    make_ready_set(root),
                    data_root=make_data_root(root),
                    narestraints_builder=lambda *_: None,
                )

    def test_completed_postmr_reports_clear_non_overwrite_message(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            source.write_text(postmr_model_text())
            run = root / "run"
            report_path = make_postmr_report(run, source)
            report = json.loads(report_path.read_text())
            report["stage"] = "postmr"
            report["status"] = "POSTMR_READY"
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(
                PostMRPreparationError, "PostMR is already complete"
            ):
                prepare_postmr(
                    run,
                    make_ready_set(root),
                    data_root=make_data_root(root),
                    narestraints_builder=lambda *_: None,
                )


if __name__ == "__main__":
    unittest.main()
