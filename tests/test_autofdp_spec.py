"""`nasolve autofdp` is specified but not implemented (docs/autofdp.md).

The command must exist with its documented options, refuse to run, and write
nothing; every engine function is a stub.
"""

import contextlib
import inspect
import io
import tempfile
import unittest
from pathlib import Path

import nasolve.autofdp as autofdp
from nasolve.cli import build_parser, main


class AutoFDPSpecificationTests(unittest.TestCase):
    def test_command_lists_the_specified_options(self):
        help_text = build_parser()._subparsers._group_actions[0].choices["autofdp"].format_help()
        for option in ("--from", "--site", "--element", "--candidates", "--theory",
                       "--tolerance", "--max-iterations", "--cycles", "--b-start",
                       "--b-step", "--occupancy-floor", "--map-sigma",
                       "--no-confirm-b", "--final", "--plan"):
            self.assertIn(option, help_text)
        self.assertIn("NOT IMPLEMENTED", help_text)

    def test_command_refuses_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "AutoMR" / "run_001"
            run.mkdir(parents=True)
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                code = main(["autofdp", str(run), "--site", "B:101", "--plan"])
            self.assertEqual(code, 1)
            self.assertIn("not implemented", stderr.getvalue())
            self.assertEqual(list(run.iterdir()), [])

    def test_engine_functions_are_documented_stubs(self):
        names = [name for name in autofdp.__all__
                 if inspect.isfunction(getattr(autofdp, name))]
        self.assertIn("execute_autofdp", names)
        for name in names:
            function = getattr(autofdp, name)
            with self.subTest(function=name):
                self.assertTrue(inspect.getdoc(function))
                arguments = [None] * len(
                    [p for p in inspect.signature(function).parameters.values()
                     if p.default is inspect.Parameter.empty
                     and p.kind is p.POSITIONAL_OR_KEYWORD]
                )
                keywords = {
                    p.name: None for p in inspect.signature(function).parameters.values()
                    if p.default is inspect.Parameter.empty and p.kind is p.KEYWORD_ONLY
                }
                with self.assertRaises(NotImplementedError):
                    function(*arguments, **keywords)

    def test_default_settings_match_the_specification(self):
        settings = autofdp.AutoFDPSettings()
        self.assertEqual(settings.tolerance, 1.0)
        self.assertEqual(settings.cycles, 3)
        self.assertEqual(settings.b_start, "surrounding")
        self.assertTrue(settings.confirm_b)
        self.assertFalse(settings.final)


if __name__ == "__main__":
    unittest.main()
