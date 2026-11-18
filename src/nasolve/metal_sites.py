"""Heavy metals and the bases around them, read from PDB coordinates.

AutoFDP finds its metal sites and their neighbouring nucleotides with this
module instead of calling NARestraints, so it has no NARestraints
dependency. The two searches are ported unchanged from NARestraints
(AcetylCholine233/NARestraints ``main`` at d1bb52d, which includes the
anomalous-scatterers branch to d0a05ec):

- ``find_heavy_atoms`` from ``restraints/anomalous.py``: HETATM atoms of the
  first model (waters excluded) whose element, from the element column or
  inferred from the residue and atom names, has atomic number >= 21;
- ``find_metal_mediated_pairs`` and ``read_atoms`` from
  ``restraints/metal_pairs.py``: for each metal ion outside nucleotides, the
  two non-neighbouring nucleotides with a base N/O/S/Se within 3.6 A (closest
  pair kept).

tests/test_metal_sites.py checks that both give the same results as
NARestraints whenever NARestraints is installed. ``scan_metal_sites`` adds
NASolve's selection (heavy metals only, not Li/Be/Na/Mg/Al/K/Ca; outside
nucleotides), as PostMR's ``scan_metal_anomalous_sites`` does with
NARestraints.

Pure Python: fixed-column PDB records, no Biopython or Phenix.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# --- ported from NARestraints restraints/anomalous.py (d1bb52d) -------------

_ELEMENT_SYMBOLS = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm",
)

#: Upper-case element symbol -> atomic number. Deuterium is treated as H.
ATOMIC_NUMBERS: dict[str, int] = {
    symbol.upper(): number for number, symbol in enumerate(_ELEMENT_SYMBOLS, start=1)
}
ATOMIC_NUMBERS["D"] = 1

#: Default threshold: scandium (Z = 21) and heavier.
DEFAULT_MIN_ATOMIC_NUMBER = 21

#: Residue names treated as water and never scanned.
WATER_RESIDUES = frozenset({"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3"})

#: The phenix.refine default strategy (Phenix 2.2.1 refinement parameters).
#: ``group_anomalous`` is not part of the default, so it must be added for the
#: anomalous groups to be refined at all.
PHENIX_DEFAULT_STRATEGY = (
    "individual_sites",
    "individual_sites_real_space",
    "individual_adp",
    "occupancies",
)

ELEMENT_COLUMN = "element column"
INFERRED_FROM_ION_RESIDUE = "inferred from ion residue name"
INFERRED_FROM_ATOM_NAME = "inferred from atom name"


@dataclass(frozen=True)
class HeavyAtomSite:
    """One heavy atom found in a HETATM record."""

    element: str          # upper-case symbol, e.g. "AG"
    atom_name: str        # stripped, e.g. "AG"
    residue_name: str     # stripped, e.g. "AG"
    chain: str            # stripped; empty when the chain column is blank
    resseq: str           # stripped residue number, e.g. "101"
    icode: str            # insertion code or ""
    altloc: str           # alternate-location indicator or ""
    occupancy: float | None
    line_number: int      # 1-based line number in the PDB file
    element_source: str   # ELEMENT_COLUMN or one of the INFERRED_* values

    @property
    def resid(self) -> str:
        return f"{self.resseq}{self.icode}"

    @property
    def label(self) -> str:
        chain = self.chain or "_"
        alt = f" alt {self.altloc}" if self.altloc else ""
        return f"{chain}:{self.resid} {self.residue_name} {self.atom_name}{alt}"

    def selection(self) -> str:
        """Phenix atom selection for this atom (all alternate locations)."""
        parts = []
        if self.chain:
            parts.append(f"chain {self.chain}")
        parts.append(f"resid {self.resid}")
        parts.append(f"name {self.atom_name}")
        return " and ".join(parts)


def element_symbol(value: str) -> str | None:
    """Return the upper-case element symbol for ``value``, or None if unknown."""
    text = "".join(ch for ch in value.strip() if ch.isalpha()).upper()
    return text if text in ATOMIC_NUMBERS else None


def parse_element_list(text: str | Iterable[str]) -> frozenset[str]:
    """Parse ``"AG,HG"`` (or an iterable of symbols) into upper-case symbols.

    Raises ValueError for anything that is not an element symbol.
    """
    if isinstance(text, str):
        items = [item for item in text.replace(";", ",").replace(" ", ",").split(",")]
    else:
        items = list(text)

    symbols: set[str] = set()
    for item in items:
        item = item.strip()
        if not item:
            continue
        symbol = item.upper()
        if symbol not in ATOMIC_NUMBERS:
            raise ValueError(f"Not an element symbol: {item!r}")
        symbols.add(symbol)

    if not symbols:
        raise ValueError("No element symbols were given.")
    return frozenset(symbols)


def infer_element(name_field: str, residue_name: str) -> tuple[str | None, str]:
    """Infer an element from the 4-character atom-name field (columns 13-16).

    Returns ``(symbol or None, how it was inferred)``.
    """
    stripped = name_field.strip()
    letters = "".join(ch for ch in stripped if ch.isalpha()).upper()

    # Ion residues such as "AG" / "AG" or "ZN" / "ZN": residue name and atom
    # name agree and are an element symbol. This also resolves "CA" calcium
    # correctly, which the column convention below would read as carbon.
    if letters and letters == residue_name.strip().upper() and letters in ATOMIC_NUMBERS:
        return letters, INFERRED_FROM_ION_RESIDUE

    # PDB convention: a one-letter element is written in column 14 (so column
    # 13 is blank or a digit); a two-letter element fills columns 13-14.
    padded = name_field.ljust(4)
    if padded[0] == " " or padded[0].isdigit():
        candidate = padded[1].upper()
    else:
        candidate = padded[0:2].upper()
    if candidate in ATOMIC_NUMBERS:
        return candidate, INFERRED_FROM_ATOM_NAME

    # Fall back to the leading letters of the stripped name.
    for length in (2, 1):
        candidate = letters[:length]
        if len(candidate) == length and candidate in ATOMIC_NUMBERS:
            return candidate, INFERRED_FROM_ATOM_NAME

    return None, INFERRED_FROM_ATOM_NAME


def _parse_occupancy(line: str) -> float | None:
    try:
        return float(line[54:60])
    except ValueError:
        return None


def find_heavy_atoms(
    pdb_filename: str | Path,
    *,
    elements: Iterable[str] | None = None,
    min_atomic_number: int = DEFAULT_MIN_ATOMIC_NUMBER,
) -> tuple[list[HeavyAtomSite], list[str]]:
    """Find heavy atoms in the HETATM records of the first model.

    If ``elements`` is given, exactly those elements are reported, whatever
    their atomic number. Otherwise every element with atomic number of at
    least ``min_atomic_number`` is reported.

    Returns ``(sites, warnings)``. Warnings report HETATM atoms whose element
    could not be determined; they are skipped rather than guessed.
    """
    wanted = parse_element_list(elements) if elements is not None else None

    sites: list[HeavyAtomSite] = []
    warnings: list[str] = []
    seen_model = False

    path = Path(pdb_filename)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            record = line[0:6]

            if record == "MODEL ":
                if seen_model:
                    break  # only the first model, as for restraint generation
                seen_model = True
                continue
            if record == "ENDMDL":
                break
            if record != "HETATM":
                continue

            residue_name = line[17:20].strip()
            if residue_name.upper() in WATER_RESIDUES:
                continue

            name_field = line[12:16]
            element = element_symbol(line[76:78]) if len(line) > 76 else None
            source = ELEMENT_COLUMN
            if element is None:
                element, source = infer_element(name_field, residue_name)

            atom_name = name_field.strip()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()

            if element is None:
                warnings.append(
                    f"Line {line_number}: could not determine the element of "
                    f"HETATM {chain or '_'}:{resseq} {residue_name} {atom_name}; skipped."
                )
                continue

            if wanted is not None:
                if element not in wanted:
                    continue
            elif ATOMIC_NUMBERS[element] < min_atomic_number:
                continue

            sites.append(
                HeavyAtomSite(
                    element=element,
                    atom_name=atom_name,
                    residue_name=residue_name,
                    chain=chain,
                    resseq=resseq,
                    icode=line[26:27].strip(),
                    altloc=line[16:17].strip(),
                    occupancy=_parse_occupancy(line),
                    line_number=line_number,
                    element_source=source,
                )
            )

    return sites, warnings


# --- ported from NARestraints restraints/metal_pairs.py (d1bb52d) -----------

#: Nonmetals and metalloids; every other element is treated as a metal.
NON_METALS = frozenset({
    "H", "D", "HE", "B", "C", "N", "O", "F", "NE", "SI", "P", "S", "CL", "AR",
    "GE", "AS", "SE", "BR", "KR", "SB", "TE", "I", "XE", "AT", "RN",
})
METALS = frozenset(symbol for symbol in ATOMIC_NUMBERS if symbol not in NON_METALS)

#: Default metals for pair detection: metals with Z >= 21.
DEFAULT_PAIR_METALS = frozenset(
    symbol for symbol in METALS if ATOMIC_NUMBERS[symbol] >= DEFAULT_MIN_ATOMIC_NUMBER
)

#: Maximum metal-to-base-atom distance (Angstrom) for a contact.
DEFAULT_MAX_DISTANCE = 3.6

#: Contacts shorter than this are reported as clashes.
CLASH_DISTANCE = 1.8

#: Elements a metal ion can bind to on a base.
LIGAND_ELEMENTS = frozenset({"N", "O", "S", "SE"})

#: Sugar and phosphate atoms, which never count as base atoms.
BACKBONE_ATOMS = frozenset(
    name.replace("'", star)
    for name in (
        "P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P",
        "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
    )
    for star in ("'", "*")
)

ResidueKey = tuple[str, str]  # (chain character as in the PDB, residue number)


@dataclass(frozen=True)
class Atom:
    record: str
    name: str
    element: str | None
    residue_name: str
    chain: str          # raw chain character (" " when blank), as Biopython uses
    resseq: str         # residue number without insertion code, e.g. "12"
    altloc: str
    xyz: tuple[float, float, float]
    line_number: int

    @property
    def residue(self) -> ResidueKey:
        return (self.chain, self.resseq)


def residue_label(key: ResidueKey, residue_name: str = "") -> str:
    chain, resseq = key
    text = f"{chain.strip() or '_'}:{resseq}"
    return f"{text} {residue_name}" if residue_name else text


@dataclass(frozen=True)
class MetalContact:
    residue: ResidueKey
    residue_name: str
    atom_name: str
    distance: float

    @property
    def label(self) -> str:
        return f"{residue_label(self.residue, self.residue_name)} {self.atom_name} {self.distance:.2f} A"


@dataclass(frozen=True)
class MetalMediatedPair:
    metal: Atom
    first: MetalContact
    second: MetalContact

    @property
    def residues(self) -> tuple[ResidueKey, ResidueKey]:
        return (self.first.residue, self.second.residue)

    @property
    def metal_label(self) -> str:
        alt = f" alt {self.metal.altloc}" if self.metal.altloc else ""
        return f"{self.metal.element} {residue_label(self.metal.residue)}{alt}"

    def describe(self) -> str:
        return (
            f"{residue_label(self.first.residue, self.first.residue_name)} <-> "
            f"{residue_label(self.second.residue, self.second.residue_name)} "
            f"via {self.metal_label} "
            f"({self.first.atom_name} {self.first.distance:.2f} A, "
            f"{self.second.atom_name} {self.second.distance:.2f} A)"
        )


def read_atoms(pdb_filename: str | Path) -> list[Atom]:
    """Read ATOM/HETATM records of the first model."""
    atoms: list[Atom] = []
    seen_model = False
    with Path(pdb_filename).open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            record = line[0:6]
            if record == "MODEL ":
                if seen_model:
                    break
                seen_model = True
                continue
            if record == "ENDMDL":
                break
            if record not in ("ATOM  ", "HETATM"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                resseq = str(int(line[22:26]))
            except ValueError:
                continue
            residue_name = line[17:20].strip()
            name_field = line[12:16]
            element = element_symbol(line[76:78]) if len(line) > 76 else None
            if element is None:
                element, _ = infer_element(name_field, residue_name)
            atoms.append(
                Atom(
                    record=record.strip(),
                    name=name_field.strip(),
                    element=element,
                    residue_name=residue_name,
                    chain=line[21:22] or " ",
                    resseq=resseq,
                    altloc=line[16:17].strip(),
                    xyz=xyz,
                    line_number=line_number,
                )
            )
    return atoms


def _adjacent(first: ResidueKey, second: ResidueKey) -> bool:
    return first[0] == second[0] and abs(int(first[1]) - int(second[1])) <= 1


class _Grid:
    """Spatial hash so each metal only checks nearby atoms."""

    def __init__(self, atoms: Iterable[Atom], cell: float):
        self.cell = cell
        self.cells: dict[tuple[int, int, int], list[Atom]] = {}
        for atom in atoms:
            self.cells.setdefault(self._key(atom.xyz), []).append(atom)

    def _key(self, xyz):
        return tuple(int(math.floor(value / self.cell)) for value in xyz)

    def near(self, xyz):
        cx, cy, cz = self._key(xyz)
        for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3):
            yield from self.cells.get((cx + dx, cy + dy, cz + dz), ())


def find_metal_mediated_pairs(
    pdb_filename: str | Path,
    *,
    max_distance: float = DEFAULT_MAX_DISTANCE,
    metals: Iterable[str] | None = None,
) -> tuple[list[MetalMediatedPair], list[str]]:
    """Return ``(pairs, messages)`` for metal-mediated base pairs in the model.

    ``metals`` limits the metal elements (default: :data:`DEFAULT_PAIR_METALS`).
    Messages explain ions that bind fewer than two bases, alternative pairs
    that were not chosen, and contacts short enough to be clashes.
    """
    if max_distance <= 0:
        raise ValueError("max_distance must be positive")
    wanted = parse_element_list(metals) if metals is not None else DEFAULT_PAIR_METALS

    atoms = read_atoms(pdb_filename)

    by_residue: dict[ResidueKey, list[Atom]] = {}
    for atom in atoms:
        by_residue.setdefault(atom.residue, []).append(atom)
    nucleotides = {
        key for key, members in by_residue.items()
        if any(atom.name in ("C1'", "C1*") for atom in members)
    }

    base_atoms = [
        atom for atom in atoms
        if atom.residue in nucleotides
        and atom.name not in BACKBONE_ATOMS
        and atom.element in LIGAND_ELEMENTS
    ]
    metal_atoms = [
        atom for atom in atoms
        if atom.record == "HETATM"
        and atom.element in wanted
        and atom.residue not in nucleotides
        and atom.residue_name.upper() not in WATER_RESIDUES
    ]

    grid = _Grid(base_atoms, max_distance)
    pairs: list[MetalMediatedPair] = []
    messages: list[str] = []
    used_pairs: set[frozenset[ResidueKey]] = set()

    for metal in metal_atoms:
        closest: dict[ResidueKey, MetalContact] = {}
        for atom in grid.near(metal.xyz):
            distance = math.dist(metal.xyz, atom.xyz)
            if distance > max_distance:
                continue
            current = closest.get(atom.residue)
            if current is None or distance < current.distance:
                closest[atom.residue] = MetalContact(
                    residue=atom.residue,
                    residue_name=atom.residue_name,
                    atom_name=atom.name,
                    distance=distance,
                )

        alt = f" alt {metal.altloc}" if metal.altloc else ""
        label = f"{metal.element} {residue_label(metal.residue)}{alt}"
        contacts = sorted(closest.values(), key=lambda contact: contact.distance)

        for contact in contacts:
            if contact.distance < CLASH_DISTANCE:
                messages.append(
                    f"{label} is only {contact.distance:.2f} A from "
                    f"{residue_label(contact.residue, contact.residue_name)} "
                    f"{contact.atom_name}; the ion position clashes and should be refined."
                )

        options = sorted(
            (
                (first, second)
                for first, second in itertools.combinations(contacts, 2)
                if not _adjacent(first.residue, second.residue)
            ),
            key=lambda option: option[0].distance + option[1].distance,
        )

        if not options:
            if contacts:
                bound = ", ".join(contact.label for contact in contacts)
                messages.append(
                    f"{label} binds base atoms of {bound}, but no two non-neighbouring "
                    f"nucleotides; no metal-mediated pair."
                )
            continue

        first, second = options[0]
        key = frozenset((first.residue, second.residue))
        for other_first, other_second in options[1:]:
            messages.append(
                f"{label}: also within {max_distance:g} A of "
                f"{residue_label(other_first.residue, other_first.residue_name)} and "
                f"{residue_label(other_second.residue, other_second.residue_name)}; "
                f"kept the closer pair."
            )
        if key in used_pairs:
            messages.append(f"{label} bridges an already detected pair; not repeated.")
            continue
        used_pairs.add(key)
        pairs.append(MetalMediatedPair(metal=metal, first=first, second=second))

    return pairs, messages



# --- NASolve's metal-site selection -----------------------------------------

#: Metals whose anomalous signal is too weak to calibrate; never metal sites.
LIGHT_METALS = frozenset({"LI", "BE", "NA", "MG", "AL", "K", "CA"})


def _nucleotide_keys(path: Path) -> set[tuple[str, str]]:
    """(chain or "_", resid) of residues with C1' and one of C2'/C3'/O4'."""
    atoms: dict[tuple[str, str], set[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line[0:6].strip().upper() not in {"ATOM", "HETATM"} or len(line) < 27:
            continue
        key = (line[21:22].strip() or "_", line[22:26].strip() + line[26:27].strip())
        atoms.setdefault(key, set()).add(line[12:16].strip())
    return {
        key for key, names in atoms.items()
        if "C1'" in names and names & {"C2'", "C3'", "O4'"}
    }


def scan_metal_sites(path: Path) -> tuple[list[dict[str, object]], list[HeavyAtomSite], list[str]]:
    """Metal ions to calibrate: heavy atoms that are metals, outside nucleotides.

    ``find_heavy_atoms`` (Z >= 21), keeping NASolve's metal elements
    (``autosol_metal.METAL_SYMBOLS``) other than :data:`LIGHT_METALS` in
    residues that are not nucleotides. Returns the
    site records (site, residue, atom name, alternate location, element and
    how it was found, coordinates), the :class:`HeavyAtomSite` objects in the
    same order, and the search's warnings. The same selection as PostMR's
    ``scan_metal_anomalous_sites``, without NARestraints.
    """
    from .autosol_metal import METAL_SYMBOLS

    heavy_metals = frozenset(METAL_SYMBOLS) - LIGHT_METALS
    sites, warnings = find_heavy_atoms(path)
    nucleotides = _nucleotide_keys(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    candidates: list[dict[str, object]] = []
    selected: list[HeavyAtomSite] = []
    for site in sites:
        chain = site.chain or "_"
        if site.element not in heavy_metals or (chain, site.resid) in nucleotides:
            continue
        line = lines[site.line_number - 1]
        selected.append(site)
        candidates.append({
            "site": f"{chain}:{site.resid}",
            "residue": site.residue_name,
            "atom_name": site.atom_name,
            "alternate": site.altloc or None,
            "element": site.element,
            "element_source": site.element_source,
            "coordinates": [round(float(line[start:end]), 3)
                            for start, end in ((30, 38), (38, 46), (46, 54))],
            "kind": "metal",
        })
    return candidates, selected, warnings


__all__ = [
    "ATOMIC_NUMBERS",
    "Atom",
    "BACKBONE_ATOMS",
    "DEFAULT_MAX_DISTANCE",
    "HeavyAtomSite",
    "LIGAND_ELEMENTS",
    "LIGHT_METALS",
    "METALS",
    "MetalContact",
    "MetalMediatedPair",
    "find_heavy_atoms",
    "find_metal_mediated_pairs",
    "read_atoms",
    "residue_label",
    "scan_metal_sites",
]
