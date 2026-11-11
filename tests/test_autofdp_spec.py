"""`nasolve autofdp` is specified in docs/autofdp.md. Step 1 (B start) is
implemented and tested in test_autofdp_bstart.py; the later steps must stay
documented stubs until they are implemented.
"""

import inspect
import unittest

import nasolve.autofdp as autofdp
from nasolve.cli import build_parser


class AutoFDPSpecificationTests(unittest.TestCase):
    def test_command_lists_the_specified_options(self):
        help_text = build_parser()._subparsers._group_actions[0].choices["autofdp"].format_help()
        for option in ("--from", "--site", "--element", "--candidates", "--theory",
                       "--tolerance", "--max-iterations", "--cycles", "--b-start",
                       "--b-step", "--occupancy-floor", "--map-sigma",
                       "--no-confirm-b", "--final", "--plan"):
            self.assertIn(option, help_text)
        self.assertIn("Only step 1 is implemented", help_text)
        self.assertIn("NOT IMPLEMENTED", help_text)

    def test_engine_functions_are_documented_stubs(self):
        names = ["plan_autofdp", "execute_autofdp", "theoretical_scattering", "fdp_round",
                 "map_value_at_site", "decide_site_update", "identity_check",
                 "confirm_b_round", "final_round"]
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
