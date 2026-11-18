"""The line-preserving PHIL scanner used to edit a manual phenix.refine .eff."""

import unittest

from nasolve.phil_text import PhilSyntaxError, edit, inside, scan, starred, unquote

GUI_STYLE = '''data_manager {
  restraint_files = "/data/restraints/PDU.cif"
  model {
    file = "/data/model with space.pdb"
    type = *x_ray neutron electron reference
  }
  default_model = "/data/model with space.pdb"
}
refinement {
  refine {
    strategy = *individual_sites *individual_sites_real_space rigid_body \\
               *individual_adp group_adp tls *occupancies *group_anomalous \\
               den
    anomalous_scatterers {
      group {
        selection = chain h resid 2   # a comment
        refine = *f_prime *f_double_prime
      }
    }
  }
  gui {
    base_output_dir = /data/Structural \\
                      Biology/run
  }
}
output {
  prefix = """PdUPdU-Hg_refine"""
  title = "a \\"quoted\\" { brace }"
  serial = 77
}
'''


class ScanTests(unittest.TestCase):
    def test_paths_values_and_round_trip(self):
        items = scan(GUI_STYLE)
        by_path = {item.path: item for item in items}
        self.assertEqual(unquote(by_path["data_manager.model.file"].value), "/data/model with space.pdb")
        self.assertEqual(
            starred(by_path["refinement.refine.strategy"].value),
            ("individual_sites", "individual_sites_real_space", "individual_adp",
             "occupancies", "group_anomalous"),
        )
        self.assertEqual(by_path["refinement.refine.anomalous_scatterers.group.selection"].value,
                         "chain h resid 2")
        self.assertEqual(by_path["refinement.gui.base_output_dir"].value,
                         "/data/Structural Biology/run")
        self.assertEqual(unquote(by_path["output.prefix"].value), "PdUPdU-Hg_refine")
        self.assertEqual(unquote(by_path["output.title"].value), 'a "quoted" { brace }')
        self.assertEqual(by_path["output.serial"].value, "77")
        self.assertIsNone(unquote("None"))
        self.assertEqual(edit(GUI_STYLE), GUI_STYLE)

    def test_removing_and_replacing_whole_items(self):
        items = scan(GUI_STYLE)
        model = next(item for item in items if item.path == "data_manager.model")
        groups = next(item for item in items if item.path == "refinement.refine.anomalous_scatterers")
        strategy = next(item for item in items if item.path == "refinement.refine.strategy")
        self.assertEqual(len(inside(items, groups)), 3)   # group scope and its two assignments
        new = GUI_STYLE[strategy.start:strategy.end].replace("*group_anomalous", "group_anomalous")
        edited = edit(GUI_STYLE, [model, *inside(items, model), groups], {strategy: new})
        self.assertNotIn("model {", edited)
        self.assertNotIn("anomalous_scatterers", edited)
        self.assertIn("*occupancies group_anomalous \\\n", edited)
        self.assertIn('  default_model = "/data/model with space.pdb"\n', edited)
        self.assertEqual([item.path for item in scan(edited)][:2], ["data_manager", "data_manager.restraint_files"])
        # The rest is untouched, character for character.
        self.assertEqual(
            edited,
            GUI_STYLE[:model.start] + GUI_STYLE[model.end:strategy.start] + new
            + GUI_STYLE[strategy.end:groups.start] + GUI_STYLE[groups.end:],
        )

    def test_syntax_errors(self):
        for text, message in (
            ("a {\n  b = 1\n", "never closed"),
            ("}\n", "unmatched"),
            ("a\n", "expected"),
            ('a = "open\n', "unterminated"),
        ):
            with self.subTest(text=text):
                with self.assertRaisesRegex(PhilSyntaxError, message):
                    scan(text)


if __name__ == "__main__":
    unittest.main()
