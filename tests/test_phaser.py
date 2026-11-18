import json
import shutil
import tempfile
import unittest
from pathlib import Path

from nasolve.automr import prepare_automr
from nasolve.phaser import (
    TFZ_FLOOR,
    TFZ_THRESHOLD,
    PhaserExecutionError,
    TFZOverrideError,
    check_tfz_override,
    execute_phaser,
    parse_best_tfz,
    tfz_override_warning,
)

from .helpers import make_dataset, make_mtz_dump, make_phaser, model_text, w_model_text


VALID = {"DA", "DC"}


class PhaserExecutionTests(unittest.TestCase):
    def _preflight(self, root: Path):
        dataset = make_dataset(root / "dataset")
        return prepare_automr(
            dataset, valid_ligand_codes=VALID,
            mtz_dump_executable=make_mtz_dump(root / "tools"),
        )

    def test_non_ascii_and_special_paths_reach_phaser_unchanged(self):
        # Regression: JSON-style quoting wrote 学... into phaser.eff, which
        # PHIL keeps literally, so Phaser could not find the reflections.
        def phil_value(eff: str, key: str) -> str:
            line = next(l.strip() for l in eff.splitlines() if l.strip().startswith(key + " = "))
            quoted = line.split(" = ", 1)[1]
            self.assertTrue(quoted.startswith('"') and quoted.endswith('"'), quoted)
            out, chars, i = [], quoted[1:-1], 0
            while i < len(chars):  # libtbx.phil: only \\ and \" are escapes
                if chars[i] == "\\" and i + 1 < len(chars) and chars[i + 1] in '\\"':
                    out.append(chars[i + 1])
                    i += 2
                else:
                    out.append(chars[i])
                    i += 1
            return "".join(out)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "学习一个" / "Structural Biology" / 'FT"Hg\\set'
            preflight = self._preflight(root)
            result = execute_phaser(
                preflight.report_path, make_phaser(root), phenix_version="2.2.1"
            )
            eff = result.parameter_path.read_text(encoding="utf-8")
            self.assertNotIn("\\u", eff)
            self.assertIn("学习一个", eff)
            report = json.loads(preflight.report_path.read_text())
            hklin = Path(phil_value(eff, "hklin"))
            self.assertEqual(hklin, Path(report["inputs"]["reflections"]))
            self.assertTrue(hklin.is_file())
            self.assertTrue(Path(phil_value(eff, "pdb")).is_file())
            self.assertEqual(Path(phil_value(eff, "output_dir")), result.parameter_path.parent.resolve())

    def test_one_ensemble_uses_checker_residue_count_and_preserves_hetatm(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root)
            result = execute_phaser(
                preflight.report_path, make_phaser(root), phenix_version="1.20.1"
            )
            self.assertEqual(result.status, "MR_SUCCESS")
            self.assertEqual(result.exit_code, 0)
            eff = result.parameter_path.read_text()
            self.assertEqual(eff.count("  ensemble {"), 1)
            self.assertNotIn("\n  model =", eff)
            self.assertIn("model_id = nasolve_model", eff)
            self.assertIn("ensembles = nasolve_model", eff)
            self.assertIn("use_hetatm = True", eff)
            self.assertIn("rmsd = 1.0", eff)
            self.assertNotIn("identity =", eff)
            self.assertIn("nres = 2", eff)
            self.assertIn("copies = 1", eff)
            self.assertTrue(result.solution_pdb.is_file())
            self.assertTrue(result.solution_mtz.is_file())
            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["execution"]["phaser"]["ensemble_count"], 1)
            self.assertEqual(report["execution"]["phaser"]["model_rmsd"], 1.0)
            self.assertEqual(report["execution"]["phaser"]["composition_nres"], 2)

    def test_preflight_inputs_rebase_after_checkout_moves(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root / "collaborator")
            source_dataset = preflight.run_directory.parent.parent
            checkout_dataset = root / "checkout" / source_dataset.name
            shutil.copytree(source_dataset, checkout_dataset)
            shutil.rmtree(root / "collaborator")
            report_path = (
                checkout_dataset / "AutoMR" / preflight.run_directory.name / "report.json"
            )

            result = execute_phaser(
                report_path,
                make_phaser(root),
                phenix_version="1.20.1",
            )

            self.assertEqual(result.status, "MR_SUCCESS")
            parameters = result.parameter_path.read_text()
            self.assertIn(str(checkout_dataset.resolve()), parameters)

    def test_relocated_preflight_rejects_changed_local_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root / "collaborator")
            source_dataset = preflight.run_directory.parent.parent
            checkout_dataset = root / "checkout" / source_dataset.name
            shutil.copytree(source_dataset, checkout_dataset)
            run = checkout_dataset / "AutoMR" / preflight.run_directory.name
            (run / "Model" / "input_model.pdb").write_text("changed\n")

            with self.assertRaisesRegex(PhaserExecutionError, "recorded checksum"):
                execute_phaser(
                    run / "report.json",
                    make_phaser(root),
                    phenix_version="1.20.1",
                )

    def test_preflight_records_and_checks_authoritative_reflections(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root)
            report = json.loads(preflight.report_path.read_text())
            self.assertRegex(report["inputs"]["reflections_sha256"], r"^[0-9a-f]{64}$")
            Path(report["inputs"]["reflections"]).write_bytes(b"changed")

            with self.assertRaisesRegex(PhaserExecutionError, "reflections"):
                execute_phaser(preflight.report_path, make_phaser(root))

    def test_tfz_review_and_failure_gates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            review = self._preflight(root / "review")
            reviewed = execute_phaser(review.report_path, make_phaser(root / "review", 7.5))
            self.assertEqual((reviewed.status, reviewed.exit_code), ("MR_REVIEW", 3))

            failed = self._preflight(root / "fail")
            failure = execute_phaser(failed.report_path, make_phaser(root / "fail", 6.9))
            self.assertEqual((failure.status, failure.exit_code), ("MR_FAILED", 4))

    def test_p1_shunt_uses_three_composition_and_search_copies(self):
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
            preflight = prepare_automr(
                dataset,
                frame_override="W",
                pair_override="D:T",
                frames_dir=root / "frames",
                allow_p1_standard=True,
                mtz_dump_executable=make_mtz_dump(root, "P 1", 1),
                valid_ligand_codes={"1AP", "DT", "DC", "DG"},
            )
            result = execute_phaser(preflight.report_path, make_phaser(root))
            eff = result.parameter_path.read_text()
            self.assertIn("num = 3", eff)
            self.assertIn("copies = 3", eff)
            self.assertEqual(eff.count("  ensemble {"), 1)

    def test_missing_tfz_or_outputs_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            missing_tfz = self._preflight(root / "tfz")
            result = execute_phaser(
                missing_tfz.report_path, make_phaser(root / "tfz", None)
            )
            self.assertEqual(result.status, "MR_FAILED")
            self.assertIn("did not report a TFZ", result.message)

            missing_outputs = self._preflight(root / "outputs")
            result = execute_phaser(
                missing_outputs.report_path,
                make_phaser(root / "outputs", 8.4, create_outputs=False),
            )
            self.assertEqual(result.status, "MR_FAILED")
            self.assertIn("both PDB and MTZ", result.message)

    def test_existing_phaser_directory_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root)
            (preflight.run_directory / "Phaser").mkdir()
            with self.assertRaisesRegex(PhaserExecutionError, "refusing to overwrite"):
                execute_phaser(preflight.report_path, make_phaser(root))

    def test_parser_prefers_best_solu_set(self):
        text = "\n".join([
            "intermediate TFZ=99.0",
            "SOLU SET RFZ=5.0 TFZ=7.2 LLG=50",
            "SOLU SET RFZ=8.0 TFZ==8.6 LLG=125.5",
        ])
        self.assertEqual(parse_best_tfz(text), (8.6, 125.5))


def low_tfz_run(root: Path, *, tfz=5.1, status="MR_FAILED", return_code=0,
                requested="F", used="F", outputs=("pdb", "mtz")) -> tuple[Path, dict]:
    """An AutoMR run as phaser.py records it, with a chosen Phaser outcome."""
    run = root / "dataset" / "AutoMR" / "run_001"
    phaser = run / "Phaser"
    phaser.mkdir(parents=True)
    record = {"return_code": return_code, "best_tfz": tfz, "best_llg": 40.0,
              "reflection_data": {"requested": requested, "labin": "FP,SIGFP",
                                  "phaser_reported_labin": "F=FP SIGF=SIGFP",
                                  "phaser_reported_kind": used},
              "solution_pdb": None, "solution_mtz": None}
    for suffix in outputs:
        path = phaser / f"mr_solution.{suffix}"
        path.write_text("fixture\n")
        record[f"solution_{suffix}"] = str(path)
    report = {"workflow": "automr", "stage": "phaser", "status": status,
              "message": f"Best TFZ {tfz} fixture", "execution": {"phaser": record}}
    (run / "report.json").write_text(json.dumps(report))
    return run, report


class TFZOverrideTests(unittest.TestCase):
    """--tfz-override accepts a Phaser result whose only problem is its TFZ."""

    def test_cut_offs_are_unchanged(self):
        self.assertEqual((TFZ_FLOOR, TFZ_THRESHOLD), (7.0, 8.0))

    def test_low_tfz_is_accepted_and_described(self):
        for tfz, status, below in ((5.1, "MR_FAILED", "below the 7.0 floor"),
                                   (7.4, "MR_REVIEW", "below the 8.0 threshold")):
            with self.subTest(tfz=tfz), tempfile.TemporaryDirectory() as directory:
                run, report = low_tfz_run(Path(directory), tfz=tfz, status=status)
                record = check_tfz_override(report, run)
                self.assertEqual(record["accepted_by"], "--tfz-override")
                self.assertEqual(record["mr_status"], status)
                self.assertEqual(record["mr_message"], report["message"])
                self.assertEqual(record["best_tfz"], tfz)
                self.assertEqual((record["tfz_floor"], record["tfz_threshold"]), (7.0, 8.0))
                self.assertIn("accepted_utc", record)
                warning = tfz_override_warning(record)
                self.assertIn(f"{status} (best TFZ {tfz:.2f}, {below})", warning)
                self.assertIn("--tfz-override", warning)

    def test_tfz_at_the_threshold_needs_no_override(self):
        with tempfile.TemporaryDirectory() as directory:
            run, report = low_tfz_run(Path(directory), tfz=8.0, status="MR_SUCCESS")
            self.assertIsNone(check_tfz_override(report, run))

    def test_outputs_are_found_by_the_run_fallback_after_a_move(self):
        with tempfile.TemporaryDirectory() as directory:
            run, report = low_tfz_run(Path(directory))
            report["execution"]["phaser"]["solution_pdb"] = "/moved/away/mr_solution.pdb"
            report["execution"]["phaser"]["solution_mtz"] = "/moved/away/mr_solution.mtz"
            self.assertIsNotNone(check_tfz_override(report, run))

    def test_any_other_failure_is_refused(self):
        cases = (
            ({"return_code": 1}, "exited with status 1"),
            ({"tfz": None}, "reported no TFZ"),
            ({"used": "I"}, "Phaser used I data although F data were requested"),
            ({"outputs": ("pdb",)}, "no solution map"),
            ({"outputs": ()}, "no solution model"),
        )
        for change, message in cases:
            with self.subTest(change=change), tempfile.TemporaryDirectory() as directory:
                run, report = low_tfz_run(Path(directory), **change)
                with self.assertRaisesRegex(TFZOverrideError, message):
                    check_tfz_override(report, run)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(TFZOverrideError, "no Phaser record"):
                check_tfz_override({"status": "MR_FAILED"}, Path(directory))


if __name__ == "__main__":
    unittest.main()
