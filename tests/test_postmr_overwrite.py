"""nasolve postmr --force-overwrite: rerun PostMR on an AutoMR run that already
has one, deleting the old PostMR and everything made from it (user's choice,
2026-11-19: "delete and restart")."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.cli import build_parser
from nasolve.postmr import PostMRPreparationError, prepare_postmr

from .test_postmr import make_data_root, make_ready_set
from .test_postmr_guessed_pairs import TRIANGLE, frameless_run

KEPT = ("Phaser/mr_solution.pdb", "CootGUI/automr/launch.json", "AutoSolMetal/run_001/report.json")
DELETED = ("PostMR", "AutoSol", "AutoRefine", "RefineDoctor",
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
    postmr(run, root)
    for relative in ("AutoRefine/checkpoints.json", "AutoRefine/round_001/refined_001.pdb",
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
            with self.assertRaisesRegex(PostMRPreparationError, "--force-overwrite deletes it"):
                postmr(run, root)
            self.assertTrue((run / "AutoRefine" / "round_001").is_dir())

    def test_deletes_downstream_reruns_and_records_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = refined_run(root)
            old_created = json.loads((run / "report.json").read_text())["postmr"]["created_utc"]
            result = postmr(run, root, force_overwrite=True)
            self.assertEqual(result.status, "POSTMR_READY")
            self.assertEqual(result.overwritten, DELETED)
            for relative in ("AutoSol", "AutoRefine", "RefineDoctor", "CootGUI/postmr",
                             "CootGUI/autosol", "CootGUI/autorefine"):
                self.assertFalse((run / relative).exists(), relative)
            for relative in KEPT:
                self.assertTrue((run / relative).is_file(), relative)
            report = json.loads((run / "report.json").read_text())
            self.assertEqual((report["stage"], report["status"]), ("postmr", "POSTMR_READY"))
            for key in ("autorefine", "autosol", "refine_doctor"):
                self.assertNotIn(key, report)
            self.assertIn("metal_merge", report)
            self.assertNotEqual(report["postmr"]["created_utc"], old_created)
            self.assertEqual(report["phaser_result"]["status"], "MR_SUCCESS")
            [record] = report["postmr_overwrites"]
            self.assertEqual(record["deleted"], list(DELETED))
            self.assertEqual((record["previous_stage"], record["previous_status"]),
                             ("autorefine", "AUTOREFINE_READY"))
            self.assertEqual(record["previous_postmr_created_utc"], old_created)
            # It can be done again; the history grows.
            postmr(run, root, force_overwrite=True)
            self.assertEqual(len(json.loads((run / "report.json").read_text())["postmr_overwrites"]), 2)

    def test_phaser_gate_is_applied_again_before_anything_is_deleted(self):
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
            self.assertTrue((run / "AutoRefine" / "round_001").is_dir())   # nothing deleted
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
            # The partial PostMR/ is refused without the flag and replaced with it.
            with self.assertRaisesRegex(PostMRPreparationError, "already exists"):
                postmr(run, root)
            self.assertEqual(postmr(run, root, force_overwrite=True).overwritten, ("PostMR",))

    def test_cli_flag(self):
        self.assertFalse(build_parser().parse_args(["postmr", "run_001"]).force_overwrite)
        self.assertTrue(build_parser().parse_args(["postmr", "run_001", "--force-overwrite"]).force_overwrite)


if __name__ == "__main__":
    unittest.main()
