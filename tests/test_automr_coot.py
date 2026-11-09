import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.cli import main
from nasolve.config import AppConfig
from nasolve.coot_runtime import CootDiscoveryError, CootInstallation
from nasolve.coot_view import launch_coot_view
from nasolve.phenix_runtime import PhenixInstallation

from .helpers import make_dataset, make_mtz_dump, make_phaser


class AutoMRCootTests(unittest.TestCase):
    def _phenix(self, root: Path, **phaser) -> PhenixInstallation:
        return PhenixInstallation(
            root=None, setup_script=None, version="2.2.1-6174",
            executables={
                "phenix.phaser": make_phaser(root, **phaser),
                "phenix.mtz.dump": make_mtz_dump(root),
            },
            environment=None, source="test",
        )

    def _main(self, argv, phenix, coot=None, launch=None, coot_error=None):
        output, errors = io.StringIO(), io.StringIO()
        discover = (
            {"side_effect": coot_error} if coot_error is not None else {"return_value": coot}
        )
        with (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli.discover_phenix", return_value=phenix),
            patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
            patch("nasolve.cli.remember_coot"),
            patch("nasolve.cli.discover_coot", **discover),
            patch("nasolve.cli.launch_coot_view", side_effect=launch) as launched,
            redirect_stdout(output), redirect_stderr(errors),
        ):
            code = main(argv)
        return code, output.getvalue(), errors.getvalue(), launched

    def test_execute_opens_phaser_model_and_map_in_coot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset")
            coot = root / "coot"
            coot.write_text("#!/bin/sh\n")
            launch = partial(launch_coot_view, launcher=lambda *a, **k: SimpleNamespace(pid=4321))
            code, out, errors, launched = self._main(
                ["automr", str(dataset), "--execute", "--show"], self._phenix(root),
                coot=CootInstallation(coot, "1.3.3", "test"), launch=launch,
            )
            self.assertEqual(code, 0, out + errors)
            self.assertIn("Opened Coot (PID 4321) with the Phaser model and map", out)
            self.assertEqual(launched.call_args.kwargs["stage"], "automr")
            run = dataset / "AutoMR" / "run_001"
            record = json.loads((run / "CootGUI" / "automr" / "launch.json").read_text())
            self.assertEqual(Path(record["model_path"]).name, "mr_solution.pdb")
            self.assertEqual(Path(record["map_path"]).name, "mr_solution.mtz")
            phaser = run / "Phaser"
            self.assertEqual(
                (phaser / "mr_solution.pdb").read_bytes(), (phaser / "PHASER.1.pdb").read_bytes()
            )

    def test_failed_score_still_opens_when_phaser_wrote_a_solution(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset")
            launch = partial(launch_coot_view, launcher=lambda *a, **k: SimpleNamespace(pid=7))
            code, out, _, launched = self._main(
                ["automr", str(dataset), "--execute", "--show"], self._phenix(root, tfz=4.4),
                coot=CootInstallation(root / "coot", "1.3.3", "test"), launch=launch,
            )
            self.assertEqual(code, 4)
            self.assertIn("MR_FAILED", out)
            launched.assert_called_once()

    def test_default_missing_outputs_and_coot_errors_do_not_open_or_change_the_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset")
            # Default: no Coot, only the command to open it later.
            code, out, _, launched = self._main(
                ["automr", str(dataset), "--execute"], self._phenix(root)
            )
            self.assertEqual(code, 0)
            self.assertIn("Open in Coot with: nasolve show", out)
            self.assertIn("--stage automr", out)
            launched.assert_not_called()

            no_outputs = root / "no_outputs"
            no_outputs.mkdir()
            code, out, _, launched = self._main(
                ["automr", str(dataset), "--execute", "--show"],
                self._phenix(no_outputs, create_outputs=False),
            )
            self.assertEqual(code, 4)
            self.assertNotIn("nasolve show", out)
            launched.assert_not_called()

            code, out, errors, launched = self._main(
                ["automr", str(dataset), "--execute", "--show"], self._phenix(root),
                coot_error=CootDiscoveryError("Coot was not found"),
            )
            self.assertEqual(code, 0)
            self.assertIn("Coot view error: Coot was not found", errors)
            self.assertIn("--stage automr", out)
            launched.assert_not_called()

    def test_one_letter_data_flags_reach_the_preflight(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset")
            for argv, expected in ((["-I"], ("I", "IMEAN,SIGIMEAN")),
                                   (["-F"], ("F", "F,SIGF")),
                                   (["--phaser-data", "I"], ("I", "IMEAN,SIGIMEAN"))):
                code, out, errors, _ = self._main(
                    ["automr", str(dataset), *argv], self._phenix(root)
                )
                self.assertEqual(code, 0, out + errors)
                run = sorted((dataset / "AutoMR").iterdir())[-1]
                record = json.loads((run / "report.json").read_text())["phaser_data"]
                self.assertEqual((record["requested"], record["labin"]), expected)
            with self.assertRaises(SystemExit), redirect_stderr(io.StringIO()):
                main(["automr", str(dataset), "-I", "-F"])

    def test_show_needs_execute_and_coot_needs_show(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset")
            code, _, errors, _ = self._main(
                ["automr", str(dataset), "--show"], self._phenix(root)
            )
            self.assertEqual(code, 2)
            self.assertIn("--show applies only with --execute", errors)
            code, _, errors, _ = self._main(
                ["automr", str(dataset), "--execute", "--coot", "/x/coot"], self._phenix(root)
            )
            self.assertEqual(code, 2)
            self.assertIn("--coot applies only with --show", errors)
            self.assertFalse((dataset / "AutoMR").exists())


if __name__ == "__main__":
    unittest.main()
