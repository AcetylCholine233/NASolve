"""DYI (O4-triazolyl-dU): the curated entry, its reviewed dictionary and placement.

DYI is dC with the N4 amino group replaced by a 1,2,4-triazol-1-yl ring (N4 is
the triazole N1). Its dictionary is a reviewed restraint dictionary used for
refinement in place of ReadySet's, which ReadySet would build "from scratch"
with every chiral centre set to "both" (DYI is not a CCD component).
"""

import json
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from nasolve.curated_ligands import (
    CURATED_LIGANDS,
    curated_dictionary,
    ligand_definition,
    validate_curated_dictionary,
)
from nasolve.ligand_profiles import AUTHORITATIVE_CODES, effective_restraints
from nasolve.checkpoints import _initial_restraints
from nasolve.postmr import (
    MutationAction,
    _restore_shared_parent_coordinates,
    build_mutation_plan,
    prepare_postmr,
)
from nasolve.residue_aliases import resolve_pair

from .helpers import make_postmr_report, pdb_record
from .test_phosphate import phosphate
from .test_postmr import PROTECTED_PARENT_ATOMS, postmr_model_text

DICTIONARY = Path(__file__).parents[1] / "src" / "nasolve" / "data" / "ligands" / "DYI.cif"
TRIAZOLE = ("N4", "N2T", "C3T", "N4T", "C5T")


def component(path: Path = DICTIONARY) -> dict:
    text = path.read_text(encoding="utf-8")
    import io
    return dict(MMCIF2Dict(io.StringIO(text[text.index("data_comp_DYI"):])))


def coordinates(values: dict) -> dict[str, tuple[float, float, float]]:
    return {
        name: (float(x), float(y), float(z))
        for name, x, y, z in zip(values["_chem_comp_atom.atom_id"], values["_chem_comp_atom.x"],
                                 values["_chem_comp_atom.y"], values["_chem_comp_atom.z"])
    }


def angle(a, b, c) -> float:
    u = [p - q for p, q in zip(a, b)]
    v = [p - q for p, q in zip(c, b)]
    dot = sum(p * q for p, q in zip(u, v))
    return math.degrees(math.acos(dot / math.dist(a, b) / math.dist(c, b)))


def dihedral(a, b, c, d) -> float:
    import numpy as np
    p0, p1, p2, p3 = (np.array(point) for point in (a, b, c, d))
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return math.degrees(math.atan2(np.cross(b1, v) @ w, v @ w))


class DYIEntryTests(unittest.TestCase):
    def test_curated_entry_and_reviewed_dictionary(self):
        ligand = ligand_definition("DYI")
        self.assertEqual(ligand.parent_code, "DC")
        self.assertIsNone(ligand.deposition_code)   # not a CCD component
        self.assertIn("DYI", AUTHORITATIVE_CODES)
        self.assertEqual(curated_dictionary("DYI"), DICTIONARY.resolve())
        validate_curated_dictionary("DYI", DICTIONARY)
        first, second = resolve_pair("DYI:G", valid_codes={"DG"})
        self.assertEqual((first.ligand_code, second.ligand_code), ("DYI", "DG"))
        values = component()
        self.assertEqual(values["_chem_comp_atom.type_energy"][
            values["_chem_comp_atom.atom_id"].index("N4")], "NR5")
        elements: dict[str, int] = {}
        for symbol in values["_chem_comp_atom.type_symbol"]:
            elements[symbol] = elements.get(symbol, 0) + 1
        # C11 H14 N5 O7 P, with the phosphate as the dianion (as in the DC dictionary).
        self.assertEqual(elements, {"C": 11, "H": 12, "N": 5, "O": 7, "P": 1})

    def test_narestraints_label_is_a_dc_analogue_narestraints_knows(self):
        try:
            from restraints.residue_library import load_residue_records
        except ImportError:
            self.skipTest("NARestraints is not installed")
        label = CURATED_LIGANDS["DYI"].narestraints_label
        record = next(r for r in load_residue_records() if str(r["Ligand code"]) == label)
        self.assertEqual((record["Base Analog"], record["Sugar Type"]), ("C", "DNA"))

    def test_ideal_coordinates_satisfy_the_restraints(self):
        values = component()
        xyz = coordinates(values)
        for a, b, target in zip(values["_chem_comp_bond.atom_id_1"], values["_chem_comp_bond.atom_id_2"],
                                values["_chem_comp_bond.value_dist"]):
            with self.subTest(bond=(a, b)):
                self.assertAlmostEqual(math.dist(xyz[a], xyz[b]), float(target), delta=0.02)
        for a, b, c, target in zip(values["_chem_comp_angle.atom_id_1"], values["_chem_comp_angle.atom_id_2"],
                                   values["_chem_comp_angle.atom_id_3"], values["_chem_comp_angle.value_angle"]):
            with self.subTest(angle=(a, b, c)):
                self.assertAlmostEqual(angle(xyz[a], xyz[b], xyz[c]), float(target), delta=2.0)
        # The triazole is coplanar with the base, N2T on the C5 side.
        self.assertAlmostEqual(abs(dihedral(*(xyz[n] for n in ("C5", "C4", "N4", "N2T")))), 0.0, delta=1.0)
        self.assertAlmostEqual(abs(dihedral(*(xyz[n] for n in ("N3", "C4", "N4", "C5T")))), 0.0, delta=1.0)
        self.assertGreater(math.dist(xyz["H5"], xyz["H5T"]), 3.5)
        # The sugar keeps dC's chirality.
        import numpy as np
        for centre, a, b, c, sign in zip(values["_chem_comp_chir.atom_id_centre"],
                                         values["_chem_comp_chir.atom_id_1"],
                                         values["_chem_comp_chir.atom_id_2"],
                                         values["_chem_comp_chir.atom_id_3"],
                                         values["_chem_comp_chir.volume_sign"]):
            if sign == "both":
                continue
            origin = np.array(xyz[centre])
            volume = np.dot(np.array(xyz[a]) - origin,
                            np.cross(np.array(xyz[b]) - origin, np.array(xyz[c]) - origin))
            self.assertEqual(volume < 0, sign == "negative", centre)

    def test_validation_rejects_an_incomplete_or_wrong_dictionary(self):
        text = DICTIONARY.read_text(encoding="utf-8")
        bad = {
            # What ReadySet makes from scratch: sugar chirality not restrained.
            "chirality": text.replace('"C1\'" "O4\'" N1    "C2\'" negative',
                                      '"C1\'" "O4\'" N1    "C2\'" both'),
            "triazole plane": "".join(line for line in text.splitlines(keepends=True)
                                      if not line.startswith("DYI plan-2 N2T ")),
            "energy types": text.replace("_chem_comp_atom.type_energy", "_chem_comp_atom.not_energy"),
            "group": text.replace(" DNA 36 24 .", " NON-POLYMER 36 24 ."),
            # dC's amino group instead of the triazole.
            "amino N4": text.replace("DYI N4    C5T    SINGLE", "DYI N4    H41    SINGLE"),
        }
        with tempfile.TemporaryDirectory() as directory:
            for name, changed in bad.items():
                with self.subTest(case=name):
                    self.assertTrue(changed != text, f"{name}: the edit did not apply")
                    path = Path(directory) / f"{name}.cif"
                    path.write_text(changed, encoding="utf-8")
                    with self.assertRaises(ValueError):
                        validate_curated_dictionary("DYI", path)

    def test_reviewed_dictionary_replaces_readysets(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            curated = root / "DYI.cif"
            curated.write_text(DICTIONARY.read_text(encoding="utf-8"), encoding="utf-8")
            generated = root / "prepared_model.ligands.cif"
            generated.write_text(
                "data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.three_letter_code\n"
                "_chem_comp.name\n_chem_comp.group\nDYI DYI Unknown non-polymer\nXYZ XYZ Unknown non-polymer\n"
                "\ndata_comp_DYI\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\nDYI N4\n"
                "\ndata_comp_XYZ\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\nXYZ C1\n",
                encoding="utf-8",
            )
            paths, dictionaries = effective_restraints([curated], generated, root)
            self.assertIn(curated, dictionaries)
            supplement = next(path for path in dictionaries if path != curated)
            self.assertNotIn("data_comp_DYI", supplement.read_text(encoding="utf-8"))
            self.assertIn("data_comp_XYZ", supplement.read_text(encoding="utf-8"))


class DYIPostMRTests(unittest.TestCase):
    def test_postmr_refines_with_the_reviewed_dictionary_not_readysets(self):
        """ReadySet builds DYI from scratch (it is not in the CCD), with the sugar
        chirality set to "both"; refinement must get the reviewed DYI.cif."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = root / "model.pdb"
            original.write_text("".join(phosphate(residue="DYI")))
            run = root / "dataset" / "AutoMR" / "run_001"
            report_path = make_postmr_report(run, original)
            report = json.loads(report_path.read_text())
            report["frame"] = None
            report["post_mr_plan"] = {"sequences": {}, "standard_pair": None, "mutations": {}}
            report_path.write_text(json.dumps(report))
            from_scratch = (DICTIONARY.read_text(encoding="utf-8")
                            .replace("negative", "both").replace("1.411", "8.888"))

            def readyset(command, cwd, **kwargs):
                (cwd / "prepared_model.updated.pdb").write_text(Path(command[1]).read_text())
                (cwd / "prepared_model.ligands.cif").write_text(from_scratch)
                return SimpleNamespace(returncode=0, stdout="mocked ReadySet")

            def no_pairs(prepared, compatibility, pair_file, output):
                pair_file.write_text("")
                return {"retained_pair_count": 0, "retained_pairs": [], "modified_sites": ["A:12"]}

            with patch("nasolve.postmr.subprocess.run", side_effect=readyset):
                result = prepare_postmr(run, root / "ready_set", modified_pairs_only=True,
                                        modified_pair_builder=no_pairs)
            postmr = json.loads(result.report_path.read_text())
            self.assertEqual(postmr["dictionary_precedence"]["authoritative_codes"], ["DYI"])
            self.assertEqual(postmr["dictionary_precedence"]["policy"],
                             "reviewed-dictionaries-over-readyset-v1")
            selected = _initial_restraints(postmr, run)
            self.assertEqual({path.name for path in selected},
                             {"DYI.cif", "linked_phosphate.cif", "linked_phosphate.phil"})
            text = "\n".join(path.read_text(encoding="utf-8") for path in selected)
            self.assertNotIn("8.888", text)
            self.assertIn("negative", text)
            # The runtime copy folds DC's alternative sugar torsions for Phenix.
            self.assertIn("_chem_comp_tor.alt_value_angle", text)
            self.assertEqual(len(postmr["linked_phosphate_profile"]["modifications"]), 1)
            compatibility = (result.postmr_directory / "Restraints" / "narestraints_input.pdb").read_text()
            self.assertIn(" C34 ", compatibility)   # NARestraints sees a dC analogue it knows
            self.assertNotIn(" DYI ", compatibility)


    def test_autofdp_gives_phenix_refine_the_reviewed_dictionary(self):
        from nasolve.autofdp import run_readyset

        from .test_autofdp_round import make_ready_set

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            start = root / "start_model.pdb"
            start.write_text("".join(phosphate(residue="DYI")))
            result = run_readyset(start, root / "ReadySet",
                                  make_ready_set(root, generated=("DYI",)))
            self.assertEqual(result.curated_codes, ("DYI",))
            self.assertEqual([path.name for path in result.dictionaries], ["DYI.cif"])
            # ReadySet's from-scratch DYI is left out; the copy folds DC's
            # alternative sugar torsions for Phenix.
            self.assertIn("_chem_comp_tor.alt_value_angle", result.dictionaries[0].read_text())


class DYIPlacementTests(unittest.TestCase):
    def test_pair_plans_an_overlay_on_dc(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "model.pdb"
            model.write_text(postmr_model_text("DC", "DG"))
            report = make_postmr_report(root / "run", model, first="DYI", second="DG", requested="DYI:G")
            actions = build_mutation_plan(json.loads(report.read_text()), model)
            self.assertEqual(
                (actions[0].before, actions[0].after, actions[0].method, actions[0].parent_code),
                ("DC", "DYI", "coot-parent-overlap", "DC"),
            )

    def test_n4_keeps_the_triazole_bond_length_after_restoration(self):
        """Coot overlays DYI on the dC parent; NASolve then restores the shared
        atoms. N4 goes on dC's C4->N4 line at DYI's 1.41 A, not to dC's 1.33 A
        position, so its bonds to the overlaid triazole stay intact."""
        xyz = coordinates(component())
        c4, n4 = xyz["C4"], xyz["N4"]
        unit = [(b - a) / math.dist(c4, n4) for a, b in zip(c4, n4)]
        dc_n4 = tuple(a + 1.325 * u for a, u in zip(c4, unit))   # dC's C4-N4
        shared = (*PROTECTED_PARENT_ATOMS, "N1", "C2", "O2", "N3", "C4", "C5", "C6")
        parent_lines, raw_lines = [], []
        serial = 0
        for name in shared:
            serial += 1
            parent_lines.append(pdb_record("ATOM", serial, name, "DC", "A", 12, element=name[0],
                                           x=xyz[name][0], y=xyz[name][1], z=xyz[name][2]))
        parent_lines.append(pdb_record("ATOM", 99, "N4", "DC", "A", 12, element="N", x=dc_n4[0],
                                       y=dc_n4[1], z=dc_n4[2], occupancy=0.8, b_factor=31.0))
        for name in (*shared, *TRIAZOLE):   # Coot's overlay: DYI on the parent
            serial += 1
            raw_lines.append(pdb_record("HETATM", serial, name, "DYI", "A", 12, element=name[0],
                                        x=xyz[name][0], y=xyz[name][1], z=xyz[name][2]))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            parent, raw, output = root / "parent.pdb", root / "raw.pdb", root / "restored.pdb"
            parent.write_text("".join(parent_lines) + "END\n")
            raw.write_text("".join(raw_lines) + "END\n")
            action = MutationAction("A:12", "DC", "DYI", "coot-parent-overlap", parent_code="DC")
            report = _restore_shared_parent_coordinates(raw, output, (action,), {"A:12": parent},
                                                        {"DYI": DICTIONARY})
            placed = {
                line[12:16].strip(): (float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                      float(line[54:60]), float(line[60:66]))
                for line in output.read_text().splitlines() if line.startswith("HETATM")
            }
        self.assertEqual(report["A:12"]["substitutions"][0]["placement"], "canonical-parent-vector")
        self.assertAlmostEqual(math.dist(placed["C4"][:3], placed["N4"][:3]), 1.411, delta=0.002)
        self.assertAlmostEqual(math.dist(placed["N4"][:3], placed["N2T"][:3]), 1.362, delta=0.01)
        self.assertAlmostEqual(math.dist(placed["N4"][:3], placed["C5T"][:3]), 1.335, delta=0.01)
        self.assertEqual(placed["N4"][3:], (0.8, 31.0))   # occupancy and B from the parent's N4
        # Restored to dC's position (0.086 A closer to C4), N4-N2T would be ~0.05 A off.
        self.assertGreater(abs(math.dist(dc_n4, xyz["N2T"]) - 1.362), 0.04)


if __name__ == "__main__":
    unittest.main()
