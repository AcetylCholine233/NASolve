"""AutoFDP step 1: copy the starting model and set each metal's B from the
nucleotides around it (NARestraints metal-pair search). Pure Python; no
Phenix program runs in this step, and no NASolve run is needed."""

import tempfile
import unittest
from pathlib import Path

from nasolve.autofdp import AutoFDPError, prepare_start_model, surrounding_b

from .helpers import pdb_record


def _residue(serial, name, chain, number, b, atoms):
    return [
        pdb_record("ATOM", serial + index, atom, name, chain, number, element=atom[0],
                   x=x, y=y, z=z, b_factor=b)
        for index, (atom, (x, y, z)) in enumerate(atoms)
    ]


def bridged_model_text(metal_b=379.99, occupancy=0.9) -> str:
    """Ag at the origin bridging A:12 DC (N3 2.3 A, B 60) and B:4 DU (O4 2.2 A,
    B 40). B:3 DA (B 10) also touches the Ag with N6 at 3.0 A but is the chain
    neighbour of B:4, so the pair is A:12-B:4. A:20 DG (B 5) is far away."""
    lines = []
    lines += _residue(1, "DC", "A", 12, 60.0, [
        ("P", (6.0, 2.0, 0.0)), ("C1'", (5.0, 1.0, 0.0)), ("O2", (3.0, 1.5, 0.0)),
        ("N3", (2.3, 0.0, 0.0)),
    ])
    lines += _residue(5, "DC", "A", 20, 5.0, [
        ("P", (20.0, 20.0, 20.0)), ("C1'", (21.0, 20.0, 20.0)),
    ])
    lines += _residue(7, "DA", "B", 3, 10.0, [
        ("P", (0.0, 6.0, 0.0)), ("C1'", (0.0, 5.0, 0.0)), ("N6", (0.0, 3.0, 0.0)),
    ])
    lines += _residue(10, "DU", "B", 4, 40.0, [
        ("P", (-6.0, -2.0, 0.0)), ("C1'", (-5.0, -1.0, 0.0)), ("O2", (-4.5, 0.5, 0.0)),
        ("O4", (-2.2, 0.0, 0.0)),
    ])
    lines.append(pdb_record("HETATM", 14, "AG", "AG", "B", 101, occupancy=occupancy,
                            element="AG", b_factor=metal_b))
    return "".join(lines) + "END\n"


def write_model(root: Path, text: str | None = None, name: str = "refined_001.pdb") -> Path:
    path = root / name
    path.write_text(text or bridged_model_text())
    return path


def _b(line: str) -> float:
    return float(line[60:66])


class BStartTests(unittest.TestCase):
    def test_metal_b_set_to_mean_of_bridged_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = write_model(root)
            original = source.read_bytes()
            output = root / "AutoFDP" / "run_001"
            result = prepare_start_model(source, output_directory=output)

            (site,) = result.sites
            self.assertEqual((site.site, site.atom_name, site.element), ("B:101", "AG", "AG"))
            self.assertEqual(site.b_before, 379.99)
            self.assertEqual(site.b_after, 50.0)   # 4 atoms at 60, 4 at 40; B:3 not used
            self.assertEqual(site.occupancy, 0.9)  # occupancy is not touched
            self.assertIn("bridged base pair", site.neighbour_b.method)
            self.assertEqual(
                [(n.residue, n.residue_name, n.mean_b) for n in site.neighbour_b.neighbours],
                [("B:4", "DU", 40.0), ("A:12", "DC", 60.0)],
            )

            # The given model is untouched, the copy byte-identical, and only
            # the metal's B column changed in the start model.
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(result.source_copy, output / "source_model.pdb")
            self.assertEqual(result.source_copy.read_bytes(), original)
            before = source.read_text().splitlines()
            after = result.start_model.read_text().splitlines()
            changed = [i for i, (x, y) in enumerate(zip(before, after)) if x != y]
            self.assertEqual(len(before), len(after))
            self.assertEqual(len(changed), 1)
            line = after[changed[0]]
            self.assertTrue(line.startswith("HETATM") and " AG " in line)
            self.assertEqual(_b(line), 50.0)
            self.assertEqual(line[:60], before[changed[0]][:60])
            self.assertEqual(line[66:], before[changed[0]][66:])

    def test_plan_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = write_model(root)
            result = prepare_start_model(source)
            self.assertEqual(result.sites[0].b_after, 50.0)
            self.assertIsNone(result.start_model)
            self.assertEqual(sorted(path.name for path in root.iterdir()), ["refined_001.pdb"])

    def test_site_selection_and_other_b_modes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = write_model(root)
            with self.assertRaisesRegex(AutoFDPError, "No metal site B:102.*metal sites: B:101"):
                prepare_start_model(source, sites=["B:102"])
            kept = prepare_start_model(source, b_start="model")
            self.assertEqual(kept.sites[0].b_after, 379.99)
            self.assertIsNone(kept.sites[0].neighbour_b)
            fixed = prepare_start_model(source, b_start="123.4", output_directory=root / "out")
            self.assertEqual(_b(next(line for line in fixed.start_model.read_text().splitlines()
                                     if line.startswith("HETATM"))), 123.4)
            for bad in ("hot", "0", "1000"):
                with self.subTest(b_start=bad), self.assertRaises(AutoFDPError):
                    prepare_start_model(source, b_start=bad)

    def test_model_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            text = "".join(l for l in bridged_model_text().splitlines(keepends=True)
                           if not l.startswith("HETATM"))
            with self.assertRaisesRegex(AutoFDPError, "No metal HETATM"):
                prepare_start_model(write_model(root, text, "no_metal.pdb"))
            with self.assertRaisesRegex(AutoFDPError, "needs a PDB-format model"):
                prepare_start_model(write_model(root, name="model.cif"))
            with self.assertRaisesRegex(AutoFDPError, "No model file"):
                prepare_start_model(root / "missing.pdb")


class NeighbourFallbackTests(unittest.TestCase):
    def _model(self, root: Path, metal_xyz) -> tuple[Path, int]:
        text = "".join(_residue(1, "DT", "A", 5, 30.0, [
            ("C1'", (5.0, 0.0, 0.0)), ("O4", (2.5, 0.0, 0.0)), ("P", (7.0, 0.0, 0.0)),
        ]) + _residue(4, "DA", "A", 9, 70.0, [
            ("C1'", (0.0, 0.0, 5.5)), ("P", (0.0, 0.0, 9.0)),
        ]))
        x, y, z = metal_xyz
        text += pdb_record("HETATM", 9, "HG", "HG", "C", 1, element="HG", x=x, y=y, z=z, b_factor=200.0)
        path = root / "model.pdb"
        path.write_text(text + "END\n")
        return path, text.count("\n")

    def test_single_base_contact_then_nearby_nucleotides_then_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path, line = self._model(root, (0.0, 0.0, 0.0))
            single = surrounding_b(path, line)   # O4 of A:5 at 2.5 A; A:9 has no base atom
            self.assertIn("base atoms within 3.6 A", single.method)
            self.assertEqual([n.residue for n in single.neighbours], ["A:5"])
            self.assertEqual(single.b_factor, 30.0)

            path, line = self._model(root, (0.0, 0.0, 3.0))   # O4 3.9 A: no base contact
            nearby = surrounding_b(path, line)
            self.assertIn("nucleotides within 6 A", nearby.method)
            self.assertEqual([n.residue for n in nearby.neighbours], ["A:9", "A:5"])
            self.assertEqual(nearby.b_factor, 46.0)   # atom-weighted: 3 x 30, 2 x 70

            path, line = self._model(root, (40.0, 40.0, 40.0))
            with self.assertRaisesRegex(AutoFDPError, "No nucleotide within 6 A"):
                surrounding_b(path, line)


if __name__ == "__main__":
    unittest.main()
