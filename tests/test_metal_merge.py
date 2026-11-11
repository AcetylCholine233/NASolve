"""In-Coot 'merge metal sites and save' command of the metal view.

Real Coot cannot run in the test environment, so the generated startup script
is executed against a small simulation of the Coot 1.x Python calls it uses
(names checked against the Coot 1.3.3 sources bundled with Phenix 2.2.1-6174).
"""

import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.coot_view import launch_coot_view, resolve_view_profile
from nasolve.model_assessment import file_sha256

from .test_autosol_metal import ENV, make_autosol, make_run, search


class FakeCoot:
    """Model molecules as lists of atoms: [chain, resno, ins, name, alt, occ, xyz]."""

    def __init__(self, files):
        self.molecules = {}
        self.dialogs = []
        self.menu_actions = {}
        for path in files:
            atoms = []
            for line in Path(path).read_text().splitlines():
                if line[:6] in ("ATOM  ", "HETATM"):
                    atoms.append([line[21].strip(), int(line[22:26]), "", line[12:16],
                                  "", float(line[54:60]), line])
            self.molecules[len(self.molecules)] = {"name": str(path), "atoms": atoms}
        self.next_imol = len(self.molecules)  # Coot never reuses molecule numbers

    # --- coot module -------------------------------------------------------
    def module(self):
        coot = types.ModuleType("coot")
        coot.graphics_n_molecules = lambda: self.next_imol
        coot.is_valid_model_molecule = lambda imol: imol in self.molecules
        coot.molecule_name = lambda imol: self.molecules[imol]["name"]
        coot.set_molecule_name = lambda imol, name: self.molecules[imol].update(name=name)
        coot.valid_labels = lambda *args: 0
        coot.info_dialog = self.dialogs.append

        def copy_molecule(imol):
            new = self.next_imol
            self.next_imol += 1
            self.molecules[new] = {
                "name": "Copy_of_" + self.molecules[imol]["name"],
                "atoms": [list(atom) for atom in self.molecules[imol]["atoms"]],
            }
            return new
        coot.copy_molecule = copy_molecule

        def residues(imol, chain):
            return sorted({(a[1], a[2]) for a in self.molecules[imol]["atoms"] if a[0] == chain})
        coot.chain_n_residues = lambda chain, imol: len(residues(imol, chain))
        coot.seqnum_from_serial_number = lambda imol, chain, n: residues(imol, chain)[n][0]
        coot.insertion_code_from_serial_number = lambda imol, chain, n: residues(imol, chain)[n][1]
        coot.residue_info_py = lambda imol, chain, resno, ins: [
            [[a[3], a[4]], [a[5], 20.0, "AG"], [0.0, 0.0, 0.0]]
            for a in self.molecules[imol]["atoms"] if (a[0], a[1], a[2]) == (chain, resno, ins)
        ]

        def set_atom_attributes_py(changes):
            for imol, chain, resno, ins, name, alt, key, value in changes:
                assert key == "occ"
                for atom in self.molecules[imol]["atoms"]:
                    if (atom[0], atom[1], atom[2], atom[3], atom[4]) == (chain, resno, ins, name, alt):
                        atom[5] = value
        coot.set_atom_attributes_py = set_atom_attributes_py

        def merge_molecules_py(sources, target):
            for source in sources:
                self.molecules[target]["atoms"].extend(
                    [list(atom) for atom in self.molecules[source]["atoms"]]
                )
            return [1, []]
        coot.merge_molecules_py = merge_molecules_py
        coot.close_molecule = lambda imol: self.molecules.pop(imol)

        def write_pdb_file(imol, path):
            lines = [atom[6][:54] + f"{atom[5]:6.2f}" + atom[6][60:]
                     for atom in self.molecules[imol]["atoms"]]
            Path(path).write_text("\n".join(lines) + "\nEND\n")
            return 0
        coot.write_pdb_file = write_pdb_file
        return coot

    def modules(self):
        utils = types.ModuleType("coot_utils")
        utils.chain_ids = lambda imol: sorted({a[0] for a in self.molecules[imol]["atoms"]})
        gui = types.ModuleType("coot_gui")
        gui.attach_module_menu_button = lambda name: name
        gui.add_simple_action_to_menu = (
            lambda menu, label, action, callback: self.menu_actions.update({(menu, label): callback})
        )
        glib = SimpleNamespace(timeout_add=lambda delay, function: function())
        repository = types.ModuleType("gi.repository")
        repository.GLib = glib
        gi = types.ModuleType("gi")
        gi.repository = repository
        return {"coot": self.module(), "coot_utils": utils, "coot_gui": gui,
                "gi": gi, "gi.repository": repository}


def _rewrite_record(run: Path, result, record) -> None:
    """Rewrite a search report and re-index its checksum in the AutoMR run."""
    result.report_path.write_text(json.dumps(record))
    report_path = run / "report.json"
    report = json.loads(report_path.read_text())
    report["metal_autosol"]["searches"][-1]["report"]["sha256"] = file_sha256(result.report_path)
    report_path.write_text(json.dumps(report))


def _site_record(root: Path, coordinates=(0.0, 0.0, 0.0), occupancy=1.5):
    run = make_run(root)
    result = search(run, make_autosol(root, [("AG", coordinates)]), root)
    sites = result.heavy_atom_model
    # AutoSol writes relative occupancies, often above 1.
    sites.write_text("\n".join(
        line[:54] + f"{occupancy:6.2f}" + line[60:] if line.startswith("HETATM") else line
        for line in sites.read_text().splitlines()
    ) + "\n")
    record = json.loads(result.report_path.read_text())
    record["outputs"]["heavy_atom_model"]["sha256"] = file_sha256(sites)
    _rewrite_record(run, result, record)
    return run, result


class MetalMergeTests(unittest.TestCase):
    def _launch(self, root: Path, target: Path):
        coot = root / "coot"
        coot.write_text("#!/bin/sh\n")
        return launch_coot_view(
            target, coot, environment=ENV,
            launcher=lambda *args, **kwargs: SimpleNamespace(pid=99),
        )

    def test_merge_command_saves_a_new_merged_model_without_overwriting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, result = _site_record(root)
            metal = result.output_directory
            view = self._launch(root, metal)
            self.assertEqual(view.metal_merge.output_directory, metal)
            self.assertEqual(view.metal_merge.merged_stem, "mr_solution-merged")
            script = view.command[view.command.index("--script") + 1]
            record = json.loads((view.working_directory / "launch.json").read_text())
            self.assertEqual(record["metal_merge"]["merged_name"], "mr_solution-merged.pdb")

            fake = FakeCoot([view.model_path, view.extra_model_paths[0]])
            namespace = {"__name__": "__main__"}
            with patch.dict(sys.modules, fake.modules()):
                exec(compile(Path(script).read_text(), script, "exec"), namespace)
                merge = fake.menu_actions[("NASolve", "Merge metal sites and save")]
                merge(None, None)
                first = metal / "mr_solution-merged.pdb"
                self.assertTrue(first.is_file())
                merged = first.read_text()
                # Starting model atoms plus the metal site, occupancy capped at 1.
                self.assertIn(" N3 ", merged)
                silver = next(line for line in merged.splitlines() if " AG " in line)
                self.assertEqual(silver[54:60], "  1.00")
                self.assertIn("occupancy 1.50 set to 1.00", fake.dialogs[-1])
                self.assertIn(str(first), fake.dialogs[-1])
                # The open originals are untouched.
                self.assertEqual(fake.molecules[1]["atoms"][0][5], 1.5)

                second = namespace["nasolve_merge_metal_sites"]()
                self.assertEqual(Path(second), metal / "mr_solution-merged-2.pdb")
                self.assertEqual(first.read_text(), merged)

    def test_sites_deleted_in_coot_and_closed_molecules(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, result = _site_record(root, occupancy=0.8)
            view = self._launch(root, result.output_directory)
            script = view.command[view.command.index("--script") + 1]
            fake = FakeCoot([view.model_path, view.extra_model_paths[0]])
            namespace = {"__name__": "__main__"}
            with patch.dict(sys.modules, fake.modules()):
                exec(compile(Path(script).read_text(), script, "exec"), namespace)
                fake.molecules[1]["atoms"].clear()  # user deleted the site
                output = Path(namespace["nasolve_merge_metal_sites"]())
                self.assertNotIn(" AG ", output.read_text())
                self.assertIn("0 metal-site atom(s)", fake.dialogs[-1])
                fake.molecules.pop(1)  # user closed the sites molecule
                self.assertIsNone(namespace["nasolve_merge_metal_sites"]())
                self.assertIn("no longer open", fake.dialogs[-1])

    def test_site_next_to_a_symmetry_copy_is_flagged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, result = _site_record(root)
            record = json.loads(result.report_path.read_text())
            record["sites"][0]["coordinates"] = [40.0, 0.0, 0.0]
            _rewrite_record(run, result, record)
            profile = resolve_view_profile(result.output_directory)
            self.assertEqual(len(profile.metal_merge.warnings), 1)
            self.assertIn("symmetry copy", profile.metal_merge.warnings[0])
            self.assertEqual(
                resolve_view_profile(run, stage="metal").metal_merge.output_directory,
                result.output_directory,
            )

    def test_other_views_have_no_merge_command(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, _ = _site_record(root)
            view = self._launch(root, run)  # default view: Phaser model and map
            self.assertIsNone(view.metal_merge)
            script = view.command[view.command.index("--script") + 1]
            self.assertNotIn("NASOLVE_MERGE", Path(script).read_text())


if __name__ == "__main__":
    unittest.main()
