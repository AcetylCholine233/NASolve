"""`nasolve autofdp` is specified in docs/autofdp.md. Step 1 (B start) is
implemented and tested in test_autofdp_bstart.py; the preflight and step 2
(one f''-only round from a model and an MTZ) in test_autofdp_round.py; step 3
(occupancy adjustment) in test_autofdp_occupancy.py; the final round
(--final-eff, partially implemented; it also serves as the B confirmation)
in test_autofdp_final.py. The identity check and the separate B-only round
must stay documented stubs until they are implemented.
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
                       "--acceptance_absolute", "--max-iterations", "--cycles", "--b-start",
                       "--b-step", "--occupancy-floor", "--map-sigma",
                       "--no-confirm-b", "--final-eff", "--plan", "--no-coot", "--coot",
                       "--phenix-root"):
            self.assertIn(option, help_text)
        flowing = " ".join(help_text.split())  # argparse wraps lines
        self.assertIn("Steps 1 to 3 are implemented", flowing)
        self.assertIn("Coot is not opened for these rounds: the maps of an f''-only round "
                      "are not diagnostic", flowing)
        self.assertIn("With --final-eff, a converged run ends with a final round", flowing)
        self.assertIn("The final round is PARTIALLY IMPLEMENTED: it also serves as the B "
                      "confirmation", flowing)
        self.assertIn("The identity check is NOT IMPLEMENTED YET", flowing)
        self.assertIn("NOT IMPLEMENTED YET", flowing)

    def test_engine_functions_are_documented_stubs(self):
        names = ["map_value_at_site", "identity_check", "confirm_b_round"]
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
        # f'' within 10% of theory (2026-11-23; it was a fixed 1 e-).
        self.assertEqual((settings.acceptance_fraction, settings.acceptance_absolute), (0.1, None))
        self.assertEqual(settings.cycles, 3)
        self.assertEqual(settings.max_iterations, 5)
        self.assertEqual(settings.occupancy_floor, 0.05)
        self.assertEqual(settings.b_start, "surrounding")
        self.assertTrue(settings.confirm_b)
        self.assertIsNone(settings.final_eff)


if __name__ == "__main__":
    unittest.main()
