"""nasolve postmr --force-overwrite / --history / --switch: rerun PostMR on an
AutoMR run that already has one. The old PostMR and everything made from it
are kept as a numbered result in PostMRHistory/ (user, 2026-11-23; from
2026-11-19 they were deleted), and the run can switch between results."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.checkpoints import inherited_paths, initialize_registry, resolve_checkpoint
from nasolve.cli import build_parser, main
from nasolve.postmr import PostMRPreparationError, prepare_postmr
from nasolve.postmr_history import PostMRHistoryError, history, switch_to
from nasolve.model_assessment import file_sha256

from .test_postmr import make_data_root, make_ready_set
from .test_postmr_guessed_pairs import TRIANGLE, frameless_run

KEPT = ("Phaser/mr_solution.pdb", "CootGUI/automr/launch.json", "AutoSolMetal/run_001/report.json")
KEPT_AS_RESULT = ("PostMR", "AutoSol", "AutoRefine", "RefineDoctor",
                  "CootGUI/postmr", "CootGUI/autosol", "CootGUI/autorefine")


def guessed(compatibility, pair_file, output):
    pair_file.write_text("")
    output.write_text("refinement.geometry_restraints.edits {}\n")
    return {"retained_pair_count": 0, "retained_pairs": []}


def readyset_ok(command, cwd, **kwargs):
    (cwd / "prepared_model.updated.pdb").write_text(Path(command[1]).read_text())
    return SimpleNamespace(returncode=0, stdout="mocked ReadySet")


def readyset_fails(command, cwd, **kwargs):
    return SimpleNamespace(returncode=1, stdout="mocked ReadySet failure")


def postmr(run: Path, root: Path, readyset=readyset_ok, **kwargs):
    data = root / "data" if (root / "data").is_dir() else make_data_root(root)
    tool = root / "phenix.ready_set" if (root / "phenix.ready_set").exists() else make_ready_set(root)
    with patch("nasolve.postmr.subprocess.run", side_effect=readyset):
        return prepare_postmr(run, tool, data_root=data, guessed_pair_builder=guessed, **kwargs)


def refined_run(root: Path) -> Path:
    """A run with PostMR done and downstream work made from it."""
    run = frameless_run(root, TRIANGLE)
    report = json.loads((run / "report.json").read_text())
    reflections = run.parent.parent / "data.mtz"
    reflections.write_bytes(b"observations")
    report.setdefault("inputs", {})["reflections"] = str(reflections)
    (run / "report.json").write_text(json.dumps(report))
    postmr(run, root)
    initialize_registry(run)   # a real checkpoint registry rooted at this PostMR
    for relative in ("AutoRefine/round_001/refined_001.pdb",
                     "AutoSol/report.json", "RefineDoctor/run_001/report.json",
                     "CootGUI/postmr/launch.json", "CootGUI/autosol/launch.json",
                     "CootGUI/autorefine/refine-001/launch.json", *KEPT):
        path = run / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            path.write_text("{}\n")
    report = json.loads((run / "report.json").read_text())
    report.update(stage="autorefine", status="AUTOREFINE_READY", message="refined",
                  autorefine={"status": "AUTOREFINE_READY"}, autosol={"status": "AUTOSOL_READY"},
                  refine_doctor={"status": "DONE"}, metal_merge={"selected": None, "models": []})
    (run / "report.json").write_text(json.dumps(report))
    return run


class ForceOverwriteTests(unittest.TestCase):
    def test_without_the_flag_a_second_postmr_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            with self.assertRaisesRegex(PostMRPreparationError, "--force-overwrite keeps it"):
                postmr(run, root)
            self.assertTrue((run / "AutoRefine" / "round_001").is_dir())

    def test_keeps_the_previous_result_reruns_and_records_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            before = {path.relative_to(run): file_sha256(path) for folder in KEPT_AS_RESULT
                      for path in (run / folder).rglob("*") if path.is_file()}
            old = json.loads((run / "report.json").read_text())
            self.assertEqual(old["postmr_result"]["id"], "postmr-001")
            result = postmr(run, root, force_overwrite=True)
            self.assertEqual(result.status, "POSTMR_READY")
            self.assertEqual((result.kept_as, result.kept), ("postmr-001", KEPT_AS_RESULT))
            kept = run / "PostMRHistory" / "postmr-001"
            # Moved unchanged, same layout; nothing deleted.
            self.assertEqual({path.relative_to(kept): file_sha256(path) for folder in KEPT_AS_RESULT
                              for path in (kept / folder).rglob("*") if path.is_file()}, before)
            for relative in ("AutoSol", "AutoRefine", "RefineDoctor", "CootGUI/postmr",
                             "CootGUI/autosol", "CootGUI/autorefine"):
                self.assertFalse((run / relative).exists(), relative)
            for relative in KEPT:
                self.assertTrue((run / relative).is_file(), relative)
            saved = json.loads((kept / "result.json").read_text())
            self.assertEqual(saved["records"]["postmr"], old["postmr"])
            self.assertEqual(saved["records"]["autorefine"], old["autorefine"])
            self.assertEqual(saved["fields"]["stage"], "autorefine")
            report = json.loads((run / "report.json").read_text())
            self.assertEqual((report["stage"], report["status"]), ("postmr", "POSTMR_READY"))
            self.assertEqual(report["postmr_result"]["id"], "postmr-002")
            for key in ("autorefine", "autosol", "refine_doctor"):
                self.assertNotIn(key, report)
            self.assertIn("metal_merge", report)
            self.assertNotEqual(report["postmr"]["created_utc"], old["postmr"]["created_utc"])
            self.assertEqual(report["phaser_result"]["status"], "MR_SUCCESS")
            [record] = report["postmr_overwrites"]
            self.assertEqual((record["kept_as"], record["moved"]), ("postmr-001", list(KEPT_AS_RESULT)))
            self.assertEqual((record["previous_stage"], record["previous_status"]),
                             ("autorefine", "AUTOREFINE_READY"))
            self.assertEqual(record["previous_postmr_created_utc"], old["postmr"]["created_utc"])
            # Again: the second result is kept too.
            third = postmr(run, root, force_overwrite=True)
            self.assertEqual(third.kept_as, "postmr-002")
            self.assertEqual(json.loads((run / "report.json").read_text())["postmr_result"]["id"], "postmr-003")
            self.assertEqual(len(json.loads((run / "report.json").read_text())["postmr_overwrites"]), 2)

    def test_switching_between_results(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            first = json.loads((run / "report.json").read_text())
            first_files = {path.relative_to(run): file_sha256(path) for folder in KEPT_AS_RESULT
                           for path in (run / folder).rglob("*") if path.is_file()}
            postmr(run, root, force_overwrite=True)   # postmr-002 current, postmr-001 kept
            second = json.loads((run / "report.json").read_text())

            rows = history(run)
            self.assertEqual([(row["id"], row["current"]) for row in rows],
                             [("postmr-001", False), ("postmr-002", True)])
            self.assertEqual((rows[0]["autorefine_rounds"], rows[0]["current_checkpoint"]), (1, "postmr"))
            self.assertEqual(rows[1]["autorefine_rounds"], 0)
            self.assertEqual(rows[0]["fields"]["stage"], "autorefine")

            kept, restored = switch_to(run, "postmr-001")
            self.assertEqual((kept["id"], restored["id"]), ("postmr-002", "postmr-001"))
            self.assertEqual({path.relative_to(run): file_sha256(path) for folder in KEPT_AS_RESULT
                              for path in (run / folder).rglob("*") if path.is_file()}, first_files)
            report = json.loads((run / "report.json").read_text())
            for key in ("postmr", "autorefine", "autosol", "refine_doctor", "stage", "status", "message"):
                self.assertEqual(report[key], first[key], key)
            self.assertEqual(report["postmr_result"]["id"], "postmr-001")
            self.assertEqual([item["to"] for item in report["postmr_switches"]], ["postmr-001"])
            self.assertFalse((run / "PostMRHistory" / "postmr-001").exists())
            self.assertTrue((run / "PostMRHistory" / "postmr-002" / "result.json").is_file())
            # The restored checkpoints resolve again, checksums included.
            _, registry = initialize_registry(run)
            paths = inherited_paths(resolve_checkpoint(registry, "postmr"), run)
            self.assertTrue(Path(paths["model"]).is_file())

            # And back.
            switch_to(run, "postmr-002")
            self.assertEqual(json.loads((run / "report.json").read_text())["postmr"], second["postmr"])
            self.assertFalse((run / "AutoRefine" / "round_001").exists())
            for bad, message in (("postmr-002", "already the current PostMR result"),
                                 ("postmr-009", "No kept PostMR result postmr-009 .*kept: postmr-001")):
                with self.subTest(bad=bad), self.assertRaisesRegex(PostMRHistoryError, message):
                    switch_to(run, bad)

            # A new PostMR after switching keeps the current one and takes a new number.
            switch_to(run, "postmr-001")
            result = postmr(run, root, force_overwrite=True)
            self.assertEqual((result.kept_as,
                              json.loads((run / "report.json").read_text())["postmr_result"]["id"]),
                             ("postmr-001", "postmr-003"))

    def test_command_line_history_and_switch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            postmr(run, root, force_overwrite=True)
            import io
            from contextlib import redirect_stderr, redirect_stdout
            from nasolve.config import AppConfig

            def run_main(arguments):
                out, err = io.StringIO(), io.StringIO()
                with (patch("nasolve.cli.load_config", return_value=AppConfig()),
                      patch("nasolve.cli._workspace_run", return_value=run),
                      redirect_stdout(out), redirect_stderr(err)):
                    code = main(arguments)
                return code, out.getvalue(), err.getvalue()

            code, out, err = run_main(["postmr", str(run), "--history"])
            self.assertEqual(code, 0, err)
            self.assertIn("  postmr-001", out)
            self.assertIn("* postmr-002", out)
            self.assertIn("AutoRefine 1 round, current postmr", out)
            self.assertIn(f"Switch with: ./nasolve postmr {run} --switch postmr-001", out)
            code, out, err = run_main(["postmr", str(run), "--switch", "postmr-001"])
            self.assertEqual(code, 0, err)
            self.assertIn("Current PostMR result: postmr-001", out)
            self.assertIn("Kept postmr-002 in", out)
            self.assertTrue((run / "AutoRefine" / "round_001").is_dir())
            code, out, err = run_main(["postmr", str(run), "--switch", "postmr-005"])
            self.assertEqual(code, 1)
            self.assertIn("No kept PostMR result postmr-005", err)

    def test_phaser_gate_is_applied_again_before_anything_is_moved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            report = json.loads((run / "report.json").read_text())
            report["phaser_result"] = {"status": "MR_FAILED",
                                       "message": "Best TFZ 6.50 is below the 7.0 acceptance floor"}
            report["execution"]["phaser"].update(return_code=0, best_tfz=6.5,
                                                 solution_mtz=str(run / "Phaser" / "mr_solution.pdb"))
            (run / "report.json").write_text(json.dumps(report))
            before = (run / "report.json").read_text()
            with self.assertRaisesRegex(PostMRPreparationError, "not MR_FAILED; --tfz-override"):
                postmr(run, root, force_overwrite=True)
            self.assertTrue((run / "AutoRefine" / "round_001").is_dir())   # nothing moved
            self.assertEqual((run / "report.json").read_text(), before)
            result = postmr(run, root, force_overwrite=True, tfz_override=True)
            self.assertEqual(result.tfz_override_record["mr_status"], "MR_FAILED")
            self.assertEqual(json.loads((run / "report.json").read_text())["phaser_result"]["status"],
                             "MR_FAILED")

    def test_older_reports_get_the_phaser_result_from_the_phaser_record(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            report = json.loads((run / "report.json").read_text())
            del report["phaser_result"]   # written before PostMR kept it
            report["execution"]["phaser"].update(
                return_code=0, best_tfz=9.1, solution_mtz=str(run / "Phaser" / "mr_solution.pdb"),
                reflection_data={"requested": "F", "labin": "F,SIGF", "phaser_reported_labin": None,
                                 "phaser_reported_kind": None})
            (run / "report.json").write_text(json.dumps(report))
            postmr(run, root, force_overwrite=True)
            self.assertEqual(json.loads((run / "report.json").read_text())["phaser_result"],
                             {"status": "MR_SUCCESS",
                              "message": "Best TFZ 9.10 passes the 8.0 AutoMR threshold"})

    def test_a_failed_rerun_leaves_the_run_at_its_phaser_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            with self.assertRaises(PostMRPreparationError):
                postmr(run, root, readyset=readyset_fails, force_overwrite=True)
            report = json.loads((run / "report.json").read_text())
            self.assertEqual((report["stage"], report["status"]), ("phaser", "MR_SUCCESS"))
            self.assertNotIn("postmr", report)
            self.assertEqual(len(report["postmr_overwrites"]), 1)
            self.assertFalse((run / "AutoRefine").exists())
            self.assertTrue((run / "PostMRHistory" / "postmr-001" / "AutoRefine" / "round_001").is_dir())
            # The partial PostMR/ is refused without the flag and replaced with it.
            with self.assertRaisesRegex(PostMRPreparationError, "already exists"):
                postmr(run, root)
            again = postmr(run, root, force_overwrite=True)   # the partial PostMR/ is kept too
            self.assertEqual((again.kept_as, again.kept), ("postmr-002", ("PostMR",)))

    def test_cli_flag(self):
        self.assertFalse(build_parser().parse_args(["postmr", "run_001"]).force_overwrite)
        self.assertTrue(build_parser().parse_args(["postmr", "run_001", "--force-overwrite"]).force_overwrite)
        self.assertEqual(build_parser().parse_args(["postmr", "run_001", "--switch", "postmr-002"]).switch,
                         "postmr-002")
        with self.assertRaises(SystemExit), patch("sys.stderr"):
            build_parser().parse_args(["postmr", "run_001", "--switch", "postmr-002", "--force-overwrite"])


if __name__ == "__main__":
    unittest.main()
