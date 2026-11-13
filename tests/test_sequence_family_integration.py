"""Standalone opt-in sequence-family freezing and PostMR integration.

Coordinates and external-tool executables below are controlled test fixtures,
not a live Phenix/Coot chemistry validation or a published structure.
"""

import copy
import json
import shutil
import sys
import tempfile
import unittest
from hashlib import sha256
from pathlib import Path
from unittest.mock import patch

from nasolve.automr import prepare_automr
from nasolve.automr_input import AutoMRInputError, format_intent, read_intent, resolve_automr_input
from nasolve.campaigns import plan_campaign
from nasolve.checkpoint_candidate import describe_checkpoint_candidate
from nasolve.candidate_recipient_comparison import compare_checkpoint_to_recipient
from nasolve.postmr import PostMRPreparationError, build_mutation_plan, prepare_postmr
from nasolve.sequence_family import load_frozen_sequence_family
from nasolve.model_compatibility import load_model_compatibility_facts
from nasolve.search_model_comparison import load_search_model_comparison
from nasolve.sequence_reference import SequenceReferenceError
from nasolve.run_context import resolve_artifact_path

from .helpers import make_dataset, make_mtz_dump, make_ready_set, pdb_record


REFERENCE = Path(__file__).resolve().parents[1] / "src/nasolve/data/sequence_references/w-metal-scaffold.json"
FORCED_W = Path(__file__).resolve().parents[1] / "MR_frames/5W6W/5W6W_noPO4.pdb"
VALID = {"DA", "DC", "DG", "DT", "1AP", "S6G", "DF", "DE", "0DA", "0DC", "0DG", "0DT"}


def original_inventory():
    # Literal inventory from the reviewed original-scaffold audit, independently
    # specified rather than manufactured by modifying the expected target.
    codes = {"A": "DA", "C": "DC", "G": "DG", "T": "DT"}
    return {
        f"{chain}:{start + i}": codes[base]
        for chain, start, sequence in (
            ("A", 1, "GAGCAGCCTGTACGGACATCA"), ("B", 1, "CCGTACA"),
            ("C", 8, "GGCTGCT"), ("D", 1, "CTGATGT"),
        ) for i, base in enumerate(sequence)
    }


def model_text(inventory=None, metal=False):
    inventory = original_inventory() if inventory is None else inventory
    lines = []
    for serial, (site, code) in enumerate(inventory.items(), 1):
        chain, resid = site.split(":")
        # One sugar atom is enough for identity/correspondence fixtures; external
        # geometry is not being claimed or tested by this synthetic coordinate set.
        lines.append(pdb_record("ATOM", serial, "C1'", code, chain, int(resid), element="C", x=serial * 5.0))
    if metal:
        lines.append(pdb_record("HETATM", 999, "MG", "MG", "Z", 999, element="MG"))
    return "".join(lines) + "END\n"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


class SequenceFamilyIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.dataset = make_dataset(self.root / "dataset")
        (self.dataset / "search.pdb").write_text(model_text(), encoding="utf-8")

    def configure(self, *, reference="w-metal-scaffold", extra="", standard=False):
        content = "[automr]\n"
        if standard:
            content += "mode = standard\nframe = W\npair = A:T\n"
        else:
            content += "mode = nonstandard\nmodel = search.pdb\n"
        if reference is not None:
            content += f"sequence_reference = {reference}\n"
        content += extra
        (self.dataset / "nasolve.txt").write_text(content, encoding="utf-8")

    def preflight(self, *, standard=False, **kwargs):
        if standard:
            catalogue = self.root / "frames/5W6W"
            catalogue.mkdir(parents=True, exist_ok=True)
            (catalogue / "A_T.pdb").write_text(model_text(), encoding="utf-8")
            (catalogue / "C_G.pdb").write_text(model_text(), encoding="utf-8")
            return prepare_automr(
                self.dataset, frames_dir=catalogue.parent,
                mtz_dump_executable=make_mtz_dump(self.root), valid_ligand_codes=VALID, **kwargs,
            )
        return prepare_automr(
            self.dataset, valid_ligand_codes=VALID,
            mtz_dump_executable=make_mtz_dump(self.root / "tools"), **kwargs,
        )

    def completed_mr(self, result):
        run = result.run_directory
        report = read_json(result.report_path)
        model = run / "Phaser/mr_solution.pdb"
        model.parent.mkdir()
        shutil.copyfile(run / "Model/input_model.pdb", model)
        report.update(stage="phaser", status="MR_SUCCESS")
        report["execution"] = {"phaser": {"solution_pdb": str(model)}}
        write_json(result.report_path, report)
        return report, model

    def actions(self, result):
        report, model = self.completed_mr(result)
        return build_mutation_plan(report, model, run_directory=result.run_directory)

    def fake_coot(self):
        executable = self.root / "coot-fixture"
        executable.write_text(
            f"#!{sys.executable}\n"
            "import os, json\nfrom pathlib import Path\n"
            "plan = json.loads(os.environ['NASOLVE_COOT_PLAN'])\n"
            "targets = {item['site']: item['after'] for item in plan}\n"
            "out = []\n"
            "for line in Path(os.environ['NASOLVE_COOT_INPUT']).read_text().splitlines(keepends=True):\n"
            "    if line.startswith(('ATOM  ', 'HETATM')):\n"
            "        site = line[21:22].strip() + ':' + line[22:27].strip()\n"
            "        if site in targets:\n"
            "            line = line[:17] + targets[site].rjust(3) + line[20:]\n"
            "    out.append(line)\n"
            "Path(os.environ['NASOLVE_COOT_OUTPUT']).write_text(''.join(out))\n"
            "print('Controlled identity-only mutation fixture')\n",
            encoding="utf-8",
        )
        executable.chmod(0o755)
        return executable

    def test_builtin_selection_round_trips_without_changing_ordinary_defaults(self):
        self.configure()
        intent = read_intent(self.dataset / "nasolve.txt")
        resolved = resolve_automr_input(self.dataset, intent, valid_ligand_codes=VALID)
        self.assertEqual(resolved.sequence_reference, REFERENCE)
        snapshot = self.dataset / "snapshot.txt"
        snapshot.write_text(format_intent(resolved), encoding="utf-8")
        self.assertEqual(read_intent(snapshot).sequence_reference, "w-metal-scaffold")
        self.configure(reference=None)
        resolved = resolve_automr_input(self.dataset, read_intent(self.dataset / "nasolve.txt"), valid_ligand_codes=VALID)
        self.assertIsNone(resolved.sequence_reference)
        self.assertNotIn("sequence_reference", format_intent(resolved))

    def test_custom_reference_path_round_trips(self):
        custom = self.dataset / "references/explicit family.json"
        custom.parent.mkdir()
        shutil.copyfile(REFERENCE, custom)
        self.configure(reference="references/explicit family.json")
        resolved = resolve_automr_input(self.dataset, read_intent(self.dataset / "nasolve.txt"), valid_ligand_codes=VALID)
        self.assertEqual(resolved.sequence_reference, custom)
        self.assertIn("sequence_reference = references/explicit family.json", format_intent(resolved))

    def test_preflight_freezes_all_targets_without_mutating_search_model(self):
        self.configure(standard=True)
        result = self.preflight(standard=True)
        run = result.run_directory
        report = read_json(result.report_path)
        self.assertEqual(result.status, "READY_POST_MR_MUTATION")
        self.assertEqual((run / "Model/input_model.pdb").read_text(), model_text())
        family = load_frozen_sequence_family(report, run)
        self.assertEqual(len(family.codes), 42)
        self.assertEqual(report["post_mr_plan"]["allow_op3_sites"], ["D:1"])
        self.assertEqual(report["sequence_family_preflight"]["differences"], [
            {"site": "A:13", "before": "DC", "after": "DT"},
            {"site": "B:3", "before": "DG", "after": "DA"},
        ])
        self.assertEqual(
            report["sequence_family_preflight"]["comparison_status"],
            "IDENTITY_DIFFERENCES",
        )
        comparison = load_search_model_comparison(report, run)
        self.assertEqual(comparison["status"], "IDENTITY_DIFFERENCES")
        self.assertEqual(comparison["correspondence"], {
            "exact_site_set": True,
            "missing_sites": [],
            "unexpected_sites": [],
        })
        self.assertEqual(comparison["model"]["chain_count"], 4)
        self.assertEqual(comparison["model"]["polymer_residue_count"], 42)
        self.assertEqual(comparison["identity"]["mismatch_count"], 2)
        self.assertEqual(comparison["identity"]["mismatches"], [
            {
                "site": "A:13",
                "model_residue_code": "DC",
                "target_residue_code": "DT",
                "target_source": "family_reference",
                "expected_postmr_correction": True,
                "route_validated": False,
            },
            {
                "site": "B:3",
                "model_residue_code": "DG",
                "target_residue_code": "DA",
                "target_source": "family_reference",
                "expected_postmr_correction": True,
                "route_validated": False,
            },
        ])
        self.assertTrue(
            (run / "Model/search_model_comparison.json").is_file()
        )
        facts = load_model_compatibility_facts(report, run)
        self.assertEqual(facts["dimensions"]["frame_identity"]["relation"], "SAME")
        self.assertEqual(facts["dimensions"]["site_set"]["relation"], "SAME")
        self.assertEqual(
            facts["dimensions"]["residue_identity"]["relation"],
            "DIFFERENT",
        )
        self.assertEqual(
            facts["dimensions"]["residue_identity"]["mismatch_count"],
            2,
        )
        self.assertEqual(
            facts["dimensions"]["construct_family"]["relation"],
            "UNKNOWN",
        )
        self.assertEqual(
            facts["evidence"]["search_model_comparison"],
            report["search_model_comparison"],
        )
        self.assertIsNone(facts["semantics"]["donor_eligibility"])
        self.assertFalse(facts["semantics"]["automatic_reuse_authorized"])
        self.assertTrue(
            (run / "Model/model_compatibility_facts.json").is_file()
        )
        self.assertIn("C:8", family.codes)
        self.assertNotIn("C:1", family.codes)
        self.assertFalse((run / "Phaser").exists())
        self.assertFalse((run / "PostMR").exists())
        for key in ("reference", "target"):
            reference = family.record[key]
            self.assertEqual(reference["anchor"], "run")
            path = resolve_artifact_path(reference, run)
            self.assertTrue(path.is_file())
            self.assertEqual(reference["sha256"], sha256(path.read_bytes()).hexdigest())

    def test_ordinary_w_retains_only_its_existing_pair_target(self):
        self.configure(reference=None, standard=True)
        result = self.preflight(standard=True)
        report, model = self.completed_mr(result)
        self.assertNotIn("sequence_family", report["post_mr_plan"])
        self.assertNotIn("sequence_family_preflight", report)
        self.assertNotIn("search_model_comparison", report)
        facts = load_model_compatibility_facts(report, result.run_directory)
        self.assertEqual(facts["dimensions"]["site_set"]["relation"], "UNKNOWN")
        self.assertEqual(
            facts["dimensions"]["residue_identity"]["relation"],
            "UNKNOWN",
        )
        self.assertFalse((result.run_directory / "Model/sequence_reference.json").exists())
        self.assertFalse(
            (result.run_directory / "Model/search_model_comparison.json").exists()
        )
        self.assertTrue(
            (result.run_directory / "Model/model_compatibility_facts.json").is_file()
        )
        self.assertEqual(len(build_mutation_plan(report, model)), 2)

    def test_reference_bytes_are_exact_and_distinct_from_content_fingerprint(self):
        data = b" \n" + REFERENCE.read_bytes() + b"\n\t"
        custom = self.dataset / "reference.json"
        custom.write_bytes(data)
        self.configure(reference="reference.json")
        result = self.preflight()
        record = read_json(result.report_path)["post_mr_plan"]["sequence_family"]
        self.assertEqual((result.run_directory / "Model/sequence_reference.json").read_bytes(), data)
        family = load_frozen_sequence_family(read_json(result.report_path), result.run_directory)
        self.assertNotEqual(record["reference"]["sha256"], family.reference.content_sha256)

    def test_full_mutation_plan_uses_existing_routes_and_noop_sites(self):
        self.configure()
        actions = self.actions(self.preflight())
        self.assertEqual(len(actions), 42)
        self.assertEqual([(a.site, a.before, a.after, a.method) for a in actions if a.method != "none"], [
            ("A:13", "DC", "DT", "coot-mutate-base"),
            ("B:3", "DG", "DA", "coot-mutate-base"),
        ])

    def test_dataset_sequences_and_site_chemistry_have_existing_precedence(self):
        self.configure(standard=True, extra=(
            "[sequences]\nA = " + "G" * 21 + "\nB = " + "G" * 7
            + "\n[mutations]\nA:12 = D\nB:4 = Q\n"
        ))
        result = self.preflight(standard=True)
        actions = {a.site: a for a in self.actions(result)}
        self.assertEqual(actions["A:12"].after, "1AP")
        self.assertEqual(actions["A:12"].method, "coot-parent-overlap")
        self.assertEqual(actions["B:4"].after, "S6G")
        self.assertEqual(actions["A:13"].after, "DG")
        family = load_frozen_sequence_family(read_json(result.report_path), result.run_directory)
        rows = {r["site"]: r for r in family.target["sites"]}
        self.assertEqual([x["residue_code"] for x in rows["A:12"]["assignments"]], ["DA", "DG", "1AP"])

    def test_chain_c_sequence_overlay_preserves_ids_8_through_14(self):
        self.configure(extra="[sequences]\nC = AAAAAAA\n")
        actions = {a.site: a.after for a in self.actions(self.preflight())}
        self.assertEqual(actions["C:8"], "DA")
        self.assertEqual(actions["C:14"], "DA")
        self.assertNotIn("C:1", actions)

    def test_matching_model_has_no_mutation_actions(self):
        self.configure()
        values = original_inventory()
        values.update({"A:13": "DT", "B:3": "DA"})
        (self.dataset / "search.pdb").write_text(model_text(values), encoding="utf-8")
        result = self.preflight()
        self.assertEqual(read_json(result.report_path)["sequence_family_preflight"]["differences"], [])
        self.assertTrue(all(a.method == "none" for a in self.actions(result)))

    def test_source_reference_can_disappear_and_relocated_run_still_uses_snapshot(self):
        custom = self.dataset / "reference.json"
        shutil.copyfile(REFERENCE, custom)
        self.configure(reference="reference.json")
        result = self.preflight()
        report, model = self.completed_mr(result)
        custom.unlink()
        relocated = self.root / "relocated dataset"
        shutil.move(str(self.dataset), relocated)
        run = relocated / "AutoMR" / result.run_directory.name
        self.assertFalse(self.dataset.exists())
        before = {str(p.relative_to(run)): p.read_bytes() for p in run.rglob("*") if p.is_file()}
        report = read_json(run / "report.json")
        family = load_frozen_sequence_family(report, run)
        self.assertEqual(len(family.codes), 42)
        actions = build_mutation_plan(report, run / "Phaser/mr_solution.pdb", run_directory=run)
        self.assertEqual(len(actions), 42)
        self.assertEqual(before, {str(p.relative_to(run)): p.read_bytes() for p in run.rglob("*") if p.is_file()})

    def test_missing_or_invalid_reference_stops_before_run_allocation(self):
        for request in ("missing.json", "", "../outside.json", str(REFERENCE)):
            with self.subTest(request=request):
                self.configure(reference=request)
                with self.assertRaises(AutoMRInputError):
                    self.preflight()
                self.assertFalse((self.dataset / "AutoMR").exists())
        custom = self.dataset / "invalid.json"
        custom.write_text('{"schema_version": true}', encoding="utf-8")
        self.configure(reference="invalid.json")
        with self.assertRaises(AutoMRInputError):
            self.preflight()
        self.assertFalse((self.dataset / "AutoMR").exists())

    def test_duplicate_reference_keys_fail_before_allocation(self):
        custom = self.dataset / "duplicate.json"
        custom.write_text(REFERENCE.read_text().replace('"kind":', '"kind": "bad", "kind":', 1), encoding="utf-8")
        self.configure(reference="duplicate.json")
        with self.assertRaisesRegex(AutoMRInputError, "Duplicate JSON key"):
            self.preflight()
        self.assertFalse((self.dataset / "AutoMR").exists())

    def test_missing_extra_or_renumbered_polymer_is_not_silently_aligned(self):
        self.configure()
        for mode in ("missing", "extra", "renumber"):
            with self.subTest(mode=mode):
                inventory = original_inventory()
                if mode == "missing":
                    del inventory["C:8"]
                elif mode == "extra":
                    inventory["Z:1"] = "DA"
                else:
                    inventory["C:1"] = inventory.pop("C:8")
                (self.dataset / "search.pdb").write_text(model_text(inventory), encoding="utf-8")
                with self.assertRaisesRegex(AutoMRInputError, "correspondence"):
                    self.preflight()
                self.assertFalse((self.dataset / "AutoMR").exists())

    def test_ambiguous_residue_identities_are_rejected(self):
        self.configure()
        data = model_text().replace("END\n", pdb_record("ATOM", 999, "N1", "DG", "A", 12, element="N") + "END\n")
        (self.dataset / "search.pdb").write_text(data, encoding="utf-8")
        with self.assertRaisesRegex(AutoMRInputError, "ambiguous"):
            self.preflight()
        self.assertFalse((self.dataset / "AutoMR").exists())

    def test_multiple_coordinate_models_are_rejected(self):
        self.configure()
        (self.dataset / "search.pdb").write_text("MODEL        1\n" + model_text() + "ENDMDL\nMODEL        2\n" + model_text(), encoding="utf-8")
        with self.assertRaisesRegex(AutoMRInputError, "one coordinate model|ambiguous"):
            self.preflight()
        self.assertFalse((self.dataset / "AutoMR").exists())

    def test_present_but_null_contract_does_not_fall_back_to_legacy(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        report["post_mr_plan"]["sequence_family"] = None
        with self.assertRaisesRegex(PostMRPreparationError, "Malformed frozen"):
            build_mutation_plan(report, model, run_directory=result.run_directory)

    def test_explicit_run_directory_is_required_for_opted_in_planner(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        with self.assertRaisesRegex(PostMRPreparationError, "explicit run_directory"):
            build_mutation_plan(report, model)

    def test_mutated_or_missing_frozen_artifacts_stop_before_postmr(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        for role in ("reference", "target"):
            ref = report["post_mr_plan"]["sequence_family"][role]
            path = result.run_directory / ref["relative_path"]
            saved = path.read_bytes()
            for mode in ("missing", "changed"):
                with self.subTest(role=role, mode=mode):
                    if mode == "missing":
                        path.unlink()
                    else:
                        path.write_bytes(saved + b" ")
                    with self.assertRaisesRegex(PostMRPreparationError, "checksum"):
                        prepare_postmr(result.run_directory, self.root / "must-not-run")
                    self.assertFalse((result.run_directory / "PostMR").exists())
                    path.write_bytes(saved)

    def test_forged_target_with_replaced_file_hash_still_fails_semantic_check(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        ref = report["post_mr_plan"]["sequence_family"]["target"]
        path = result.run_directory / ref["relative_path"]
        target = read_json(path)
        target["sites"][0]["residue_code"] = "DC"
        write_json(path, target)
        ref.update(sha256=sha256(path.read_bytes()).hexdigest(), size=path.stat().st_size)
        with self.assertRaisesRegex(PostMRPreparationError, "differs from its reference/overlays"):
            build_mutation_plan(report, model, run_directory=result.run_directory)

    def test_changed_dataset_intent_or_mirror_flag_cannot_recompile_old_run(self):
        self.configure(standard=True)
        result = self.preflight(standard=True)
        original, model = self.completed_mr(result)
        for mode in ("sequence", "pair", "mutation", "mirror", "frame"):
            with self.subTest(mode=mode):
                report = copy.deepcopy(original)
                if mode == "sequence":
                    report["post_mr_plan"]["sequences"]["B"] = "AAAAAAA"
                elif mode == "pair":
                    report["post_mr_plan"]["standard_pair"]["ligand_codes"] = ["DC", "DG"]
                elif mode == "mutation":
                    report["post_mr_plan"]["mutations"]["A:12"] = {"requested": "D", "ligand_code": "1AP"}
                elif mode == "mirror":
                    report["inputs"]["mirror"] = True
                else:
                    report["frame"]["name"] = "3GBI"
                with self.assertRaisesRegex(PostMRPreparationError, "input intent has changed"):
                    build_mutation_plan(report, model, run_directory=result.run_directory)

    def test_bad_contract_reference_metadata_cannot_bypass_integrity(self):
        self.configure()
        result = self.preflight()
        original, model = self.completed_mr(result)
        for key, value in (("anchor", "absolute"), ("sha256", ""), ("size", True), ("relative_path", "../../escape.json")):
            with self.subTest(key=key):
                report = copy.deepcopy(original)
                report["post_mr_plan"]["sequence_family"]["reference"][key] = value
                with self.assertRaises(PostMRPreparationError):
                    build_mutation_plan(report, model, run_directory=result.run_directory)

    def test_changed_legacy_sequence_file_is_not_new_reference_authority(self):
        self.configure(reference=None, standard=True)
        result = self.preflight(standard=True)
        (result.run_directory / "Model/seq_base.txt").write_text("AAAA\n", encoding="utf-8")
        report, model = self.completed_mr(result)
        self.assertEqual(len(build_mutation_plan(report, model)), 2)
        self.assertIsNone(load_frozen_sequence_family(report, result.run_directory))

    def test_mirrored_mutation_safety_gate_remains_in_force(self):
        self.configure(extra="mirror = true\n")
        def fake_mirror(source, destination):
            lines = []
            for line in source.read_text().splitlines(keepends=True):
                if line.startswith("ATOM  "):
                    line = line[:17] + ("0" + line[17:20].strip()).rjust(3) + line[20:]
                lines.append(line)
            destination.write_text("".join(lines), encoding="utf-8")
            return destination
        result = self.preflight(mirror_transformer=fake_mirror)
        report, model = self.completed_mr(result)
        with self.assertRaisesRegex(PostMRPreparationError, "guarded.*unmirror"):
            build_mutation_plan(report, model, run_directory=result.run_directory)
        with self.assertRaisesRegex(PostMRPreparationError, "guarded.*unmirror"):
            prepare_postmr(result.run_directory, self.root / "must-not-run")
        self.assertFalse((result.run_directory / "PostMR").exists())

    def test_reference_polymer_cannot_silently_change_sugar_type(self):
        self.configure()
        data = model_text().replace("END\n", pdb_record("ATOM", 998, "O2'", "DA", "A", 12, element="O") + "END\n")
        (self.dataset / "search.pdb").write_text(data, encoding="utf-8")
        result = self.preflight()
        report, model = self.completed_mr(result)
        with self.assertRaisesRegex(PostMRPreparationError, "polymer type disagrees"):
            build_mutation_plan(report, model, run_directory=result.run_directory)

    def test_controlled_postmr_pipeline_executes_two_changes_and_audits_42_sites(self):
        self.configure()
        (self.dataset / "search.pdb").write_text(model_text(metal=True), encoding="utf-8")
        result = self.preflight()
        report, model = self.completed_mr(result)
        raw_before = model.read_bytes()
        readyset = make_ready_set(self.root)
        with patch("nasolve.sequence_family.known_ligand_codes", return_value=VALID):
            prepared = prepare_postmr(result.run_directory, readyset, coot_executable=self.fake_coot())
        payload = read_json(prepared.report_path)
        audit = payload["sequence_family_audit"]
        self.assertEqual(audit["status"], "PASS")
        self.assertEqual(audit["target_count"], 42)
        self.assertEqual(audit["model"]["anchor"], "run")
        self.assertEqual(model.read_bytes(), raw_before)
        self.assertIn(" MG ", prepared.model_path.read_text())
        changed = [a for a in payload["mutation_actions"] if a["method"] != "none"]
        self.assertEqual([a["site"] for a in changed], ["A:13", "B:3"])
        self.assertEqual(len(payload["mutation_actions"]), 42)
        self.assertEqual(read_json(result.report_path)["status"], "POSTMR_READY")

        self.assertFalse((result.run_directory / "AutoRefine").exists())
        with patch(
            "nasolve.checkpoint_candidate.known_ligand_codes",
            return_value=VALID,
        ):
            candidate = describe_checkpoint_candidate(result.run_directory)
        self.assertFalse((result.run_directory / "AutoRefine").exists())
        comparison = candidate["run_context"]["checkpoint_model_target_comparison"]
        self.assertEqual(comparison["status"], "EXACT")
        self.assertEqual(comparison["identity"]["mismatch_count"], 0)
        self.assertEqual(comparison["model"]["polymer_residue_count"], 42)
        self.assertEqual(
            candidate["run_context"]["sequence_family"]["reference"]["id"],
            "w-metal-scaffold",
        )
        self.assertEqual(
            candidate["run_context"]["source_model_provenance_verification"],
            "frozen-model-compatibility-facts",
        )
        self.assertTrue(candidate["source"]["locally_reusable"])
        self.assertIsNone(candidate["semantics"]["donor_eligibility"])
        self.assertFalse(candidate["semantics"]["automatic_reuse_authorized"])

        with patch(
            "nasolve.checkpoint_candidate.known_ligand_codes",
            return_value=VALID,
        ):
            cross = compare_checkpoint_to_recipient(
                result.run_directory,
                result.run_directory,
            )
        self.assertEqual(cross["dimensions"]["site_set"]["relation"], "SAME")
        self.assertEqual(
            cross["dimensions"]["residue_identity"]["relation"],
            "SAME",
        )
        self.assertEqual(
            cross["dimensions"]["target_reference_context"]["relation"],
            "SAME",
        )
        self.assertEqual(
            cross["dimensions"]["residue_identity"]["mismatch_count"],
            0,
        )
        self.assertIsNone(cross["semantics"]["donor_eligibility"])
        self.assertIsNone(cross["semantics"]["recipient_compatibility"])
        self.assertFalse(cross["semantics"]["rescue_authorized"])

    def test_forced_original_scaffold_normalizes_sequence_and_terminal_phosphate(self):
        self.configure(standard=True, extra="model = 5W6W_noPO4.pdb\n")
        catalogue = self.root / "frames/5W6W"
        catalogue.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(FORCED_W, catalogue / "5W6W_noPO4.pdb")
        result = self.preflight(standard=True)
        report, model = self.completed_mr(result)
        self.assertEqual(report["inputs"]["model_provider"]["kind"], "explicit-standard-model")
        self.assertEqual(report["post_mr_plan"]["allow_op3_sites"], ["D:1"])

        readyset = make_ready_set(self.root)

        def restraint_builder(model, pairs, output):
            output.write_text("geometry_restraints.edits {}\n")

        with patch("nasolve.sequence_family.known_ligand_codes", return_value=VALID):
            prepared = prepare_postmr(
                result.run_directory,
                readyset,
                coot_executable=self.fake_coot(),
                narestraints_builder=restraint_builder,
            )
        payload = read_json(prepared.report_path)
        self.assertEqual(payload["sequence_family_audit"]["status"], "PASS")
        self.assertEqual(payload["sequence_family_audit"]["target_count"], 42)

        changed = [a for a in payload["mutation_actions"] if a["method"] != "none"]
        self.assertEqual([a["site"] for a in changed], ["A:13", "B:3"])

        d1_atoms = {
            line[12:16].strip()
            for line in prepared.model_path.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
            and line[21:22].strip() == "D"
            and line[22:26].strip() == "1"
        }
        self.assertTrue({"P", "OP1", "OP2", "OP3"} <= d1_atoms)
        self.assertEqual(read_json(result.report_path)["status"], "POSTMR_READY")

    def test_prepared_model_extra_polymer_fails_complete_target_audit(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        readyset = make_ready_set(self.root)
        script = readyset.read_text()
        extra = pdb_record("ATOM", 999, "C1'", "DA", "Z", 1, element="C")
        # Append an extra synthetic polymer residue; all 42 target sites still match.
        readyset.write_text(script + "printf '%s\\n' " + __import__('shlex').quote(extra.rstrip()) + " >> prepared_model.updated.pdb\n", encoding="utf-8")
        # ReadySet's earlier atom-count guard must still reject this first.
        with self.assertRaisesRegex(PostMRPreparationError, "atom count"):
            prepare_postmr(result.run_directory, readyset, coot_executable=self.fake_coot())
        self.assertEqual(read_json(result.report_path)["status"], "MR_SUCCESS")

    def test_campaign_snapshots_selected_reference_without_silent_drop(self):
        self.configure(standard=True)
        frames = self.root / "frames/5W6W"
        frames.mkdir(parents=True)
        (frames / "C_G.pdb").write_text(model_text(), encoding="utf-8")
        result = plan_campaign(self.root, datasets=("dataset",), frames_directory=frames.parent)
        entry = result["datasets"][0]
        self.assertEqual(entry["status"], "DISCOVERED")
        self.assertEqual(entry["effective_config"]["sequence_reference"], "w-metal-scaffold")
        frozen = entry["inputs"]["sequence_reference"]
        self.assertEqual(frozen["anchor"], "campaign")
        self.assertTrue(frozen["relative_path"].startswith("NASolveCampaign/resources/"))
        self.assertEqual(
            (self.root / frozen["relative_path"]).read_bytes(),
            REFERENCE.read_bytes(),
        )
        self.assertFalse((self.dataset / "AutoMR").exists())

    def test_missing_contract_cannot_turn_an_opted_in_run_into_legacy(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        del report["post_mr_plan"]["sequence_family"]
        with self.assertRaisesRegex(PostMRPreparationError, "no frozen target contract"):
            build_mutation_plan(report, model, run_directory=result.run_directory)

    def test_postmr_revalidates_actual_mr_site_inventory_before_creating_output(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        model.write_text(model.read_text() + pdb_record("ATOM", 999, "C1'", "DA", "Z", 1, element="C"), encoding="utf-8")
        with self.assertRaisesRegex(PostMRPreparationError, "correspondence"):
            prepare_postmr(result.run_directory, self.root / "must-not-run")
        self.assertFalse((result.run_directory / "PostMR").exists())

    def test_input_and_plan_must_name_the_same_reference_artifact(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        report["inputs"]["sequence_reference"]["sha256"] = "a" * 64
        with self.assertRaisesRegex(PostMRPreparationError, "provenance disagrees"):
            build_mutation_plan(report, model, run_directory=result.run_directory)

    def test_same_atom_count_with_wrong_readyset_identity_fails_final_audit(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        from nasolve.postmr import _run_readyset
        def changed_identity(*args, **kwargs):
            value = _run_readyset(*args, **kwargs)
            output = value[0]
            lines = []
            for line in output.read_text().splitlines(keepends=True):
                if line.startswith("ATOM  ") and line[21:27].strip() == "A  13":
                    line = line[:17] + " DC" + line[20:]
                lines.append(line)
            output.write_text("".join(lines), encoding="utf-8")
            return value
        with (
            patch("nasolve.postmr._run_readyset", side_effect=changed_identity),
            self.assertRaisesRegex(PostMRPreparationError, "identity audit failed.*A:13"),
        ):
            prepare_postmr(result.run_directory, make_ready_set(self.root), coot_executable=self.fake_coot())
        self.assertEqual(read_json(result.report_path)["status"], "MR_SUCCESS")

    def test_target_artifact_mutated_during_tool_execution_is_not_accepted(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        from nasolve.postmr import _run_readyset
        def changed_contract(*args, **kwargs):
            value = _run_readyset(*args, **kwargs)
            with (result.run_directory / "Model/sequence_family_target.json").open("a") as handle:
                handle.write(" ")
            return value
        with (
            patch("nasolve.postmr._run_readyset", side_effect=changed_contract),
            self.assertRaisesRegex(PostMRPreparationError, "identity audit failed.*checksum"),
        ):
            prepare_postmr(result.run_directory, make_ready_set(self.root), coot_executable=self.fake_coot())
        self.assertEqual(read_json(result.report_path)["status"], "MR_SUCCESS")

    def test_custom_rna_reference_uses_existing_rna_mutation_route(self):
        custom = self.dataset / "rna.json"
        write_json(custom, {
            "schema_version": 1, "kind": "sequence-defined-reference", "id": "rna-example", "version": "1",
            "chains": [{"chain": "X", "residue_ids": ["10", "11"], "polymer": "RNA", "sequence": "AG"}],
        })
        self.configure(reference="rna.json")
        data = "".join(pdb_record("ATOM", i, atom, code, "X", resid, element=atom[0]) for i, (atom, code, resid) in enumerate([
            ("C1'", "A", 10), ("O2'", "A", 10), ("C1'", "C", 11), ("O2'", "C", 11),
        ], 1)) + "END\n"
        (self.dataset / "search.pdb").write_text(data, encoding="utf-8")
        actions = self.actions(self.preflight())
        self.assertEqual([(a.site, a.after, a.method) for a in actions], [
            ("X:10", "A", "none"), ("X:11", "G", "coot-mutate-base"),
        ])

    def test_new_reference_version_affects_only_new_runs(self):
        custom = self.dataset / "reference.json"
        shutil.copyfile(REFERENCE, custom)
        self.configure(reference="reference.json")
        first = self.preflight()
        first_report, first_model = self.completed_mr(first)
        updated = read_json(custom)
        updated["version"] = "experimental-2"
        updated["chains"][0]["sequence"] = "GAGCAGCCTGTACGGACATCA"
        write_json(custom, updated)
        second = self.preflight()
        second_report, second_model = self.completed_mr(second)
        first_codes = {a.site: a.after for a in build_mutation_plan(first_report, first_model, run_directory=first.run_directory)}
        second_codes = {a.site: a.after for a in build_mutation_plan(second_report, second_model, run_directory=second.run_directory)}
        self.assertEqual(first_codes["A:13"], "DT")
        self.assertEqual(second_codes["A:13"], "DC")
        self.assertNotEqual(first.run_directory, second.run_directory)

    def test_missing_native_phosphate_atoms_are_not_excused_by_sequence_coverage(self):
        self.configure(standard=True)
        result = self.preflight(standard=True)
        report, model = self.completed_mr(result)
        with self.assertRaisesRegex(PostMRPreparationError, "5'-terminal phosphate construction failed"):
            prepare_postmr(result.run_directory, self.root / "must-not-run", coot_executable=self.fake_coot())
        self.assertEqual(read_json(result.report_path)["status"], "MR_SUCCESS")



    def test_overlay_cannot_introduce_site_chemistry_absent_from_frozen_input(self):
        self.configure()
        result = self.preflight()
        report, model = self.completed_mr(result)
        contract = report["post_mr_plan"]["sequence_family"]
        contract["overlays"]["dataset_site_codes"]["A:13"] = "DG"
        with self.assertRaisesRegex(PostMRPreparationError, "site overlays disagree"):
            build_mutation_plan(report, model, run_directory=result.run_directory)


    def test_duplicate_atom_identities_do_not_pass_as_unique_correspondence(self):
        self.configure()
        data = model_text()
        duplicate = data.splitlines(keepends=True)[0]
        (self.dataset / "search.pdb").write_text(data + duplicate, encoding="utf-8")
        with self.assertRaisesRegex(AutoMRInputError, "Duplicate coordinate atom identities"):
            self.preflight()
        self.assertFalse((self.dataset / "AutoMR").exists())


if __name__ == "__main__":
    unittest.main()
