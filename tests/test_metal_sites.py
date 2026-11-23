"""nasolve.metal_sites: the heavy-atom and metal-pair searches AutoFDP uses,
ported from NARestraints so AutoFDP does not need NARestraints.

- AutoFDP runs end to end with NARestraints made unimportable.
- Whenever NARestraints is installed, the ported searches give exactly its
  results, and ``scan_metal_sites`` exactly PostMR's selection.
"""

import importlib.abc
import importlib.util
import sys
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

from nasolve import metal_sites
from nasolve.autofdp import STATUS_CONVERGED, prepare_start_model, surrounding_b

from .helpers import pdb_record
from .test_autofdp_bstart import bridged_model_text
from .test_autofdp_round import _execute, _plan, make_inputs

try:
    HAVE_NARESTRAINTS = importlib.util.find_spec("restraints") is not None
except ImportError:   # blocked, or a broken installation
    HAVE_NARESTRAINTS = False


class _RefuseNARestraints(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name == "restraints" or name.startswith("restraints."):
            raise ImportError(f"{name} is blocked for this test")
        return None


class without_narestraints:
    """Make ``import restraints`` (NARestraints) fail inside the block."""

    def __enter__(self):
        self.saved = {name: module for name, module in sys.modules.items()
                      if name == "restraints" or name.startswith("restraints.")}
        for name in self.saved:
            del sys.modules[name]
        self.finder = _RefuseNARestraints()
        sys.meta_path.insert(0, self.finder)
        return self

    def __exit__(self, *exc):
        sys.meta_path.remove(self.finder)
        sys.modules.update(self.saved)
        return False


def tricky_model_text() -> str:
    """Heavy-atom edge cases: element column blank (ion residue, atom-name
    convention), calcium written CA/CA, a water, an alternate location, a
    metal inside a nucleotide (5-iodo base), two Hg near one base pair, and a
    second MODEL that must be ignored."""
    def blank_element(line: str) -> str:
        return line[:76] + "  \n"

    lines = ["MODEL        1\n", bridged_model_text().replace("END\n", "")]
    lines += [
        blank_element(pdb_record("HETATM", 20, "HG", "HG", "C", 1, x=0.5, y=0.4, z=0.2, element="HG")),
        pdb_record("HETATM", 21, "HG", "HG", "C", 2, x=-0.4, y=0.3, z=0.1, element="HG",
                   occupancy=0.5),
        blank_element(pdb_record("HETATM", 22, "CA", "CA", "D", 1, x=30.0, y=0.0, z=0.0, element="CA")),
        pdb_record("HETATM", 23, "O", "HOH", "W", 1, x=1.0, y=1.0, z=1.0, element="O"),
        pdb_record("HETATM", 24, "ZN", "ZN", "E", 1, x=40.0, y=40.0, z=40.0, element="ZN"),
        pdb_record("HETATM", 25, "C1'", "5IU", "F", 7, x=50.0, y=0.0, z=0.0, element="C"),
        pdb_record("HETATM", 26, "C2'", "5IU", "F", 7, x=51.0, y=0.0, z=0.0, element="C"),
        pdb_record("HETATM", 27, "I5", "5IU", "F", 7, x=52.0, y=0.0, z=0.0, element="I"),
        "ENDMDL\n",
        "MODEL        2\n",
        pdb_record("HETATM", 30, "AU", "AU", "G", 1, element="AU"),
        "ENDMDL\nEND\n",
    ]
    out = []
    for line in "".join(lines).splitlines(keepends=True):
        if line.startswith("HETATM   21"):
            line = line[:16] + "A" + line[17:]   # alternate location A (column 17)
        out.append(line)
    return "".join(out)


class WithoutNARestraintsTests(unittest.TestCase):
    def test_autofdp_runs_without_narestraints(self):
        with tempfile.TemporaryDirectory() as directory, without_narestraints():
            root = Path(directory)
            with self.assertRaises(ImportError):
                import restraints  # noqa: F401
            model, mtz = make_inputs(root)
            start = prepare_start_model(model)
            self.assertEqual(start.sites[0].b_after, 50.0)
            neighbour = surrounding_b(model, start.sites[0].line_number)
            self.assertEqual(neighbour.method, "bridged base pair (metal-pair search)")
            self.assertEqual([item.residue for item in neighbour.neighbours], ["B:4", "A:12"])
            # The fixture's f'' 2.6 vs theory 2.03 converges in a fixed 1 e-
            # window (the default 10% would send it to review).
            result = _execute(_plan(root, model, mtz), root, acceptance_absolute=1.0)
            self.assertEqual(result.status, STATUS_CONVERGED, result.message)


@unittest.skipUnless(HAVE_NARESTRAINTS, "NARestraints is not installed")
class SameAsNARestraintsTests(unittest.TestCase):
    def _models(self, root: Path) -> list[Path]:
        paths = []
        for name, text in (("bridged.pdb", bridged_model_text()), ("tricky.pdb", tricky_model_text())):
            path = root / name
            path.write_text(text)
            paths.append(path)
        return paths

    def test_heavy_atoms_and_pairs_match_narestraints(self):
        from restraints import anomalous, metal_pairs

        with tempfile.TemporaryDirectory() as directory:
            for path in self._models(Path(directory)):
                with self.subTest(model=path.name):
                    ours, our_warnings = metal_sites.find_heavy_atoms(path)
                    theirs, their_warnings = anomalous.find_heavy_atoms(path)
                    self.assertEqual([asdict(site) for site in ours], [asdict(site) for site in theirs])
                    self.assertEqual(our_warnings, their_warnings)
                    self.assertEqual([asdict(atom) for atom in metal_sites.read_atoms(path)],
                                     [asdict(atom) for atom in metal_pairs.read_atoms(path)])
                    our_pairs, our_messages = metal_sites.find_metal_mediated_pairs(path)
                    their_pairs, their_messages = metal_pairs.find_metal_mediated_pairs(path)
                    self.assertEqual([pair.describe() for pair in our_pairs],
                                     [pair.describe() for pair in their_pairs])
                    self.assertEqual(our_messages, their_messages)
            self.assertTrue(our_pairs)   # the tricky model has bridged pairs to compare
            for name in ("ATOMIC_NUMBERS", "WATER_RESIDUES", "DEFAULT_MIN_ATOMIC_NUMBER"):
                self.assertEqual(getattr(metal_sites, name), getattr(anomalous, name))
            for name in ("METALS", "DEFAULT_MAX_DISTANCE", "CLASH_DISTANCE",
                         "LIGAND_ELEMENTS", "BACKBONE_ATOMS"):
                self.assertEqual(getattr(metal_sites, name), getattr(metal_pairs, name))

    def test_metal_site_selection_matches_postmr(self):
        from nasolve.postmr import scan_metal_anomalous_sites

        with tempfile.TemporaryDirectory() as directory:
            for path in self._models(Path(directory)):
                with self.subTest(model=path.name):
                    ours, our_sites, our_warnings = metal_sites.scan_metal_sites(path)
                    theirs, their_sites, their_warnings = scan_metal_anomalous_sites(path)
                    self.assertEqual(ours, theirs)
                    self.assertEqual([asdict(site) for site in our_sites],
                                     [asdict(site) for site in their_sites])
                    self.assertEqual(our_warnings, their_warnings)
            self.assertEqual([candidate["site"] for candidate in ours],
                             ["B:101", "C:1", "C:2", "E:1"])   # no CA, water, 5IU iodine, model 2


if __name__ == "__main__":
    unittest.main()
