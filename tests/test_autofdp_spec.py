"""`nasolve autofdp` is specified in docs/autofdp.md. Step 1 (B start) is
implemented and tested in test_autofdp_bstart.py; the preflight and step 2
(one f''-only round from a model and an MTZ) in test_autofdp_round.py; step 3
(occupancy adjustment) in test_autofdp_occupancy.py. The later steps must stay
documented stubs until they are implemented.
"""

import inspect
import unittest

import nasolve.autofdp as autofdp
from nasolve.cli import build_parser


class AutoFDPSpecificationTests(unittest.TestCase):
    def test_command_lists_the_specified_options(self):
        help_text = build_parser()._subparsers._group_actions[0].choices["autofdp"].format_help()
        for option in ("MODEL.pdb", "DATA.mtz", "--output", "--wavelength", "--cif",
                       "--site", "--element", "--candidates", "--theory",
                       "--tolerance", "--max-iterations", "--cycles", "--b-start",
                       "--b-step", "--occupancy-floor", "--map-sigma",
                       "--no-confirm-b", "--final", "--plan", "--phenix-root"):
            self.assertIn(option, help_text)
        self.assertNotIn("--coot", help_text)        # AutoFDP never opens Coot
        self.assertNotIn("--no-coot", help_text)     # still accepted, not shown
        flowing = " ".join(help_text.split())  # argparse wraps lines
        self.assertIn("Steps 1 to 3 are implemented", flowing)
        self.assertIn("Coot is not opened: the maps of an f''-only round are not diagnostic",
                      flowing)
        self.assertIn("NOT IMPLEMENTED YET", flowing)

    def test_engine_functions_are_documented_stubs(self):
        names = ["map_value_at_site", "identity_check", "confirm_b_round", "final_round"]
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
        self.assertEqual(settings.max_iterations, 5)
        self.assertEqual(settings.occupancy_floor, 0.05)
        self.assertEqual(settings.b_start, "surrounding")
        self.assertTrue(settings.confirm_b)
        self.assertFalse(settings.final)


if __name__ == "__main__":
    unittest.main()
