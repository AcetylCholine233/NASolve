"""Curated and compatibility residues used during guarded model preparation.

These entries are deliberately small and explicit.  NASolve does not run
eLBOW for every residue: a reviewed dictionary is bundled once and copied
into each run that needs it.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from math import dist, isfinite
from pathlib import Path


@dataclass(frozen=True)
class AtomSubstitution:
    parent_atom: str
    target_atom: str
    anchor_atom: str


@dataclass(frozen=True)
class RingSubstituent:
    target_atom: str
    anchor_atom: str
    ring_neighbors: tuple[str, str]
    minimum_bond_length: float = 1.8
    maximum_bond_length: float = 2.5


@dataclass(frozen=True)
class CuratedLigand:
    code: str
    dictionary_filename: str
    accepted_model_labels: tuple[str, ...]
    narestraints_label: str
    description: str
    parent_code: str | None = None
    deposition_code: str | None = None
    required_bonds: tuple[tuple[str, str], ...] = ()
    forbidden_bonds: tuple[tuple[str, str], ...] = ()
    atom_substitutions: tuple[AtomSubstitution, ...] = ()
    ring_substituents: tuple[RingSubstituent, ...] = ()
    # A reviewed restraint dictionary is used for refinement as it is, in
    # place of the one ReadySet generates (``ligand_profiles.AUTHORITATIVE_CODES``).
    # Each group in ``required_planes`` must lie in one of its planes.
    reviewed_restraints: bool = False
    required_planes: tuple[tuple[str, ...], ...] = ()
    # Bonds about which a reviewed dictionary must not restrain rotation: no
    # torsion about the bond, and no plane with atoms on both sides of it.
    free_rotation: tuple[tuple[str, str], ...] = ()


CURATED_LIGANDS: dict[str, CuratedLigand] = {
    "DE": CuratedLigand(
        code="DE",
        dictionary_filename="DE.cif",
        accepted_model_labels=("DE",),
        narestraints_label="DE",
        description="4-thiothymidine under the laboratory PDB-compatible label DE",
        parent_code="DT",
        required_bonds=(("C4", "S4"),),
        forbidden_bonds=(("N3", "S4"),),
        atom_substitutions=(AtomSubstitution("O4", "S4", "C4"),),
    ),
    "DF": CuratedLigand(
        code="DF",
        dictionary_filename="DF.cif",
        accepted_model_labels=("DF",),
        narestraints_label="DF",
        description="2-thiothymidine under the laboratory PDB-compatible label DF",
        parent_code="DT",
        deposition_code="A1AAZ",
        required_bonds=(("C2", "S1"),),
        forbidden_bonds=(("N3", "S1"),),
        atom_substitutions=(AtomSubstitution("O2", "S1", "C2"),),
    ),
    "1AP": CuratedLigand(
        code="1AP",
        dictionary_filename="1AP.cif",
        accepted_model_labels=("1AP",),
        narestraints_label="1AP",
        description="2,6-diaminopurine nucleotide",
        parent_code="DA",
        deposition_code="1AP",
        required_bonds=(("C1'", "N9"), ("C2", "N2"), ("C6", "N6"), ("P", "O5'")),
        forbidden_bonds=(("O4'", "N1"),),
        reviewed_restraints=True,
        required_planes=(("N1", "C2", "N2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"),),
    ),
    "S6G": CuratedLigand(
        code="S6G",
        dictionary_filename="S6G.cif",
        accepted_model_labels=("S6G",),
        narestraints_label="S6G",
        description="6-thio-2'-deoxyguanosine-5'-monophosphate",
        parent_code="DG",
        deposition_code="S6G",
        required_bonds=(("C6", "S6"),),
        forbidden_bonds=(("N1", "S6"),),
        atom_substitutions=(AtomSubstitution("O6", "S6", "C6"),),
    ),
    "C38": CuratedLigand(
        code="C38",
        dictionary_filename="C38.cif",
        accepted_model_labels=("C38",),
        narestraints_label="C38",
        description="5-iodo-2'-deoxycytidine-5'-monophosphate",
        parent_code="DC",
        deposition_code="C38",
        required_bonds=(("C5", "I"),),
        forbidden_bonds=(("C4", "I"), ("C6", "I")),
        ring_substituents=(RingSubstituent("I", "C5", ("C4", "C6")),),
    ),
    "DU": CuratedLigand(
        code="DU",
        dictionary_filename="DU.cif",
        accepted_model_labels=("DU",),
        narestraints_label="DU",
        description=(
            "2'-deoxyuridine-5'-monophosphate (official wwPDB CCD DU, as bundled "
            "with Phenix 2.2.1-6174)"
        ),
        # dU is dT without the 5-methyl: every DU atom is a DT atom, so Coot
        # overlays DU on DT and all coordinates come back from the parent.
        parent_code="DT",
        deposition_code="DU",
        required_bonds=(("C1'", "N1"), ("C2", "O2"), ("C4", "O4"), ("C5", "C6"), ("P", "O5'")),
        forbidden_bonds=(("C5", "C7"),),
    ),
    "5IU": CuratedLigand(
        code="5IU",
        dictionary_filename="5IU.cif",
        accepted_model_labels=("5IU",),
        narestraints_label="5IU",
        description="5-iodo-2'-deoxyuridine-5'-monophosphate",
        parent_code="DT",
        deposition_code="5IU",
        required_bonds=(("C5", "I5"),),
        forbidden_bonds=(("C4", "I5"), ("C6", "I5")),
        ring_substituents=(RingSubstituent("I5", "C5", ("C4", "C6")),),
    ),
    "DYI": CuratedLigand(
        code="DYI",
        dictionary_filename="DYI.cif",
        accepted_model_labels=("DYI",),
        # NARestraints has no DYI; C34 (N4-methyl-dC) is an N4-substituted dC
        # with the same base-atom names, so DYI is still a modified residue.
        narestraints_label="C34",
        description=(
            "O4-triazolyl-dU (supplier name): 4-(1,2,4-triazol-1-yl)pyrimidin-2-one "
            "2'-deoxynucleotide, under the laboratory label DYI (not a CCD component)"
        ),
        # dC bonding (N3=C4, N on C4); the triazole N1 is DYI's N4.
        parent_code="DC",
        required_bonds=(
            ("C1'", "N1"), ("C4", "N4"), ("N4", "N2T"), ("N2T", "C3T"),
            ("C3T", "N4T"), ("N4T", "C5T"), ("C5T", "N4"), ("P", "O5'"),
        ),
        forbidden_bonds=(("N4", "H41"), ("N4", "H42"), ("C4", "O4"), ("C4", "N2T")),
        # N4 stays on dC's C4->N4 line, at DYI's C4-N4 length (1.41 A, not dC's
        # 1.33 A), so the triazole Coot overlays keeps its bond to N4.
        atom_substitutions=(AtomSubstitution("N4", "N4", "C4"),),
        reviewed_restraints=True,
        # Each ring flat, with its bond to the other ring in its plane; the
        # central bond C4-N4 is free to rotate (user, 2026-11-23).
        required_planes=(
            ("N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"),
            ("C4", "N4", "N2T", "C3T", "N4T", "C5T"),
        ),
        free_rotation=(("C4", "N4"),),
    ),
}

CURATED_LIGAND_CODES = frozenset(CURATED_LIGANDS)


_PARENT_CODES = {
    ("DNA", "A"): "DA",
    ("DNA", "C"): "DC",
    ("DNA", "G"): "DG",
    ("DNA", "T"): "DT",
    ("DNA", "U"): "DU",
    ("RNA", "A"): "A",
    ("RNA", "C"): "C",
    ("RNA", "G"): "G",
    ("RNA", "U"): "U",
}


def ligand_definition(code: str) -> CuratedLigand:
    """Return a curated override or infer a conservative CCD definition."""
    if code in CURATED_LIGANDS:
        return CURATED_LIGANDS[code]
    try:
        from restraints.residue_library import load_residue_records
    except ImportError as exc:  # pragma: no cover - installation error path
        raise ValueError("NARestraints is required to infer modified residues") from exc
    matches = [
        record
        for record in load_residue_records()
        if str(record.get("Ligand code")) == code
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one NARestraints record for {code}, found {len(matches)}"
        )
    record = matches[0]
    sugar = str(record.get("Sugar Type") or "").upper()
    base = str(record.get("Base Analog") or "").upper()
    try:
        parent_code = _PARENT_CODES[(sugar, base)]
    except KeyError as exc:
        raise ValueError(
            f"Cannot infer a canonical parent for {code}: "
            f"Sugar Type={sugar!r}, Base Analog={base!r}"
        ) from exc
    return CuratedLigand(
        code=code,
        dictionary_filename=f"{code}.cif",
        accepted_model_labels=(code,),
        narestraints_label=code,
        description=str(record.get("Name") or f"CCD component {code}"),
        parent_code=parent_code,
        deposition_code=code,
    )


def ligand_data_directory(data_root: Path | None = None) -> Path:
    if data_root is not None:
        return Path(data_root) / "ligands"
    return Path(__file__).resolve().parent / "data" / "ligands"


def curated_dictionary(code: str, data_root: Path | None = None) -> Path:
    try:
        ligand = CURATED_LIGANDS[code]
    except KeyError as exc:
        raise KeyError(f"No curated NASolve ligand is registered for {code}") from exc
    path = ligand_data_directory(data_root) / ligand.dictionary_filename
    if not path.is_file():
        raise FileNotFoundError(
            f"Curated dictionary for {code} is missing: {path}"
        )
    return path


def ligand_dictionary(code: str, data_root: Path | None = None) -> Path:
    """Resolve the reviewed override or a local official CCD dictionary."""
    ligand = ligand_definition(code)
    path = ligand_data_directory(data_root) / ligand.dictionary_filename
    if not path.is_file():
        raise FileNotFoundError(
            f"No local CCD dictionary for {code}: expected {path}"
        )
    return path


def validate_ligand_dictionary(code: str, path: Path) -> None:
    """Validate curated topology or, generically, the dictionary identity."""
    if code in CURATED_LIGANDS:
        validate_curated_dictionary(code, path)
        return
    component_ids: set[str] = set()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        try:
            fields = shlex.split(raw_line, comments=False)
        except ValueError:
            continue
        if len(fields) >= 2 and fields[0] == "_chem_comp.id":
            component_ids.add(fields[1])
    if code not in component_ids:
        found = ", ".join(sorted(component_ids)) or "none"
        raise ValueError(
            f"Dictionary {path} does not declare _chem_comp.id {code}; found {found}"
        )


def _dictionary_bonds(path: Path) -> set[frozenset[str]]:
    """Read the atom pairs from a CIF ``_chem_comp_bond`` loop."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers: list[str] = []
        cursor = index + 1
        while cursor < len(lines) and lines[cursor].lstrip().startswith("_"):
            headers.append(lines[cursor].strip())
            cursor += 1
        atom_1 = "_chem_comp_bond.atom_id_1"
        atom_2 = "_chem_comp_bond.atom_id_2"
        if atom_1 not in headers or atom_2 not in headers:
            continue
        first_index = headers.index(atom_1)
        second_index = headers.index(atom_2)
        bonds: set[frozenset[str]] = set()
        while cursor < len(lines):
            row = lines[cursor].strip()
            if not row or row == "#" or row == "loop_" or row.startswith("data_"):
                break
            if row.startswith("_"):
                break
            # ``posix=False`` preserves unquoted prime atom names such as
            # O5', which occur in older eLBOW dictionaries.
            fields = shlex.split(row, comments=True, posix=False)
            if len(fields) > max(first_index, second_index):
                atoms = []
                for field_index in (first_index, second_index):
                    atom = fields[field_index]
                    if (
                        len(atom) >= 2
                        and atom[0] == atom[-1]
                        and atom[0] in {"'", '"'}
                    ):
                        atom = atom[1:-1]
                    atoms.append(atom)
                bonds.add(frozenset(atoms))
            cursor += 1
        return bonds
    raise ValueError(f"Dictionary has no _chem_comp_bond loop: {path}")


def dictionary_atom_names(path: Path) -> dict[str, str]:
    """Atom name -> element from a CIF ``_chem_comp_atom`` loop ({} if none)."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers: list[str] = []
        cursor = index + 1
        while cursor < len(lines) and lines[cursor].lstrip().startswith("_"):
            headers.append(lines[cursor].strip())
            cursor += 1
        if "_chem_comp_atom.atom_id" not in headers:
            continue
        name_index = headers.index("_chem_comp_atom.atom_id")
        element_index = (headers.index("_chem_comp_atom.type_symbol")
                         if "_chem_comp_atom.type_symbol" in headers else None)
        names: dict[str, str] = {}
        while cursor < len(lines):
            row = lines[cursor].strip()
            if not row or row in {"#", "loop_"} or row.startswith(("data_", "_")):
                break
            fields = shlex.split(row, comments=True, posix=False)
            if len(fields) > name_index:
                name = fields[name_index]
                if len(name) >= 2 and name[0] == name[-1] and name[0] in {"'", '"'}:
                    name = name[1:-1]
                element = (fields[element_index].upper()
                           if element_index is not None and len(fields) > element_index else "")
                names[name] = element
            cursor += 1
        return names
    return {}


def dictionary_ideal_bond_length(path: Path, atom_1: str, atom_2: str) -> float:
    """Return an atom-pair distance from the preferred ideal CIF coordinates."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers: list[str] = []
        cursor = index + 1
        while cursor < len(lines) and lines[cursor].lstrip().startswith("_"):
            headers.append(lines[cursor].strip())
            cursor += 1
        atom_header = "_chem_comp_atom.atom_id"
        coordinate_header_sets = (
            (
                "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
                "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
                "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
            ),
            (
                "_chem_comp_atom.x",
                "_chem_comp_atom.y",
                "_chem_comp_atom.z",
            ),
        )
        if atom_header not in headers:
            continue
        coordinate_headers = next(
            (
                coordinate_set
                for coordinate_set in coordinate_header_sets
                if all(header in headers for header in coordinate_set)
            ),
            None,
        )
        if coordinate_headers is None:
            continue
        atom_index = headers.index(atom_header)
        coordinate_indices = tuple(headers.index(header) for header in coordinate_headers)
        coordinates: dict[str, tuple[float, float, float]] = {}
        while cursor < len(lines):
            row = lines[cursor].strip()
            if not row or row == "#" or row == "loop_" or row.startswith("data_"):
                break
            if row.startswith("_"):
                break
            fields = shlex.split(row, comments=True, posix=False)
            if len(fields) > max(atom_index, *coordinate_indices):
                atom = fields[atom_index].strip("'\"")
                if atom in {atom_1, atom_2}:
                    try:
                        coordinates[atom] = tuple(
                            float(fields[field_index])
                            for field_index in coordinate_indices
                        )
                    except ValueError as exc:
                        raise ValueError(
                            f"Dictionary has invalid ideal coordinates for {atom}: {path}"
                        ) from exc
            cursor += 1
        if atom_1 in coordinates and atom_2 in coordinates:
            return dist(coordinates[atom_1], coordinates[atom_2])
        raise ValueError(
            f"Dictionary ideal-coordinate loop is missing {atom_1} or {atom_2}: {path}"
        )
    raise ValueError(f"Dictionary has no usable ideal atom coordinates: {path}")


#: Sugar chiral centres a reviewed DNA dictionary must restrain with a sign.
_SUGAR_CHIRAL_CENTRES = ("C1'", "C3'", "C4'")


def _validate_reviewed_restraints(ligand: CuratedLigand, path: Path) -> None:
    """A reviewed dictionary replaces ReadySet's, so it must be complete.

    The official CCD graph alone is not a refinement restraint dictionary.
    """
    from .ligand_profiles import _blocks

    code = ligand.code
    values: dict[str, list[str]] = {}
    for block in _blocks(path):
        for key, value in block.values.items():
            if key == "data_":
                continue
            if key in values:
                raise ValueError(f"{code} dictionary has repeated category field {key}")
            values[key] = value
    if values.get("_chem_comp.id") != [code] or any(
        set(values.get(key, [])) != {code}
        for key in ("_chem_comp_atom.comp_id", "_chem_comp_bond.comp_id")
    ):
        raise ValueError(f"{code} dictionary has conflicting component identities")
    names = values.get("_chem_comp_atom.atom_id", [])
    energies = values.get("_chem_comp_atom.type_energy", [])
    if not names or len(energies) != len(names) or any(v in ("", ".", "?") for v in energies):
        raise ValueError(f"{code} requires complete nonbonded energy types, not a raw CCD graph")
    if values.get("_chem_comp.group") != ["DNA"]:
        raise ValueError(f"{code} restraint dictionary must use the reviewed DNA component group")
    for key in ("_chem_comp_bond.value_dist", "_chem_comp_bond.value_dist_esd"):
        raw = values.get(key, [])
        try:
            valid = len(raw) == len(values["_chem_comp_bond.atom_id_1"]) and bool(raw) and all(
                isfinite(float(v)) and float(v) > 0 for v in raw)
        except (ValueError, KeyError):
            valid = False
        if not valid:
            raise ValueError(f"{code} requires finite positive numerical restraints: {key}")
    planes: dict[str, set[str]] = {}
    for plane_id, atom in zip(values.get("_chem_comp_plane_atom.plane_id", []),
                              values.get("_chem_comp_plane_atom.atom_id", [])):
        planes.setdefault(plane_id, set()).add(atom)
    for required in ligand.required_planes:
        if not any(set(required) <= members for members in planes.values()):
            raise ValueError(f"{code} plane definition is incomplete: {', '.join(required)}")
    if ligand.free_rotation:
        neighbours: dict[str, set[str]] = {}
        for first, second in zip(values.get("_chem_comp_bond.atom_id_1", []),
                                 values.get("_chem_comp_bond.atom_id_2", [])):
            neighbours.setdefault(first, set()).add(second)
            neighbours.setdefault(second, set()).add(first)

        def side(start: str, blocked: str) -> set[str]:
            found, todo = {start}, [start]
            while todo:
                for atom in neighbours.get(todo.pop(), ()):
                    if atom != blocked and atom not in found:
                        found.add(atom)
                        todo.append(atom)
            return found

        torsions = list(zip(values.get("_chem_comp_tor.atom_id_2", []),
                            values.get("_chem_comp_tor.atom_id_3", [])))
        for first, second in ligand.free_rotation:
            if any({b, c} == {first, second} for b, c in torsions):
                raise ValueError(
                    f"{code} must not restrain rotation about {first}-{second}: torsion about it")
            one, other = side(first, second) - {first}, side(second, first) - {second}
            for plane_id, members in planes.items():
                if members & one and members & other:
                    raise ValueError(
                        f"{code} must not restrain rotation about {first}-{second}: "
                        f"plane {plane_id} spans both sides")
    signs = dict(zip(values.get("_chem_comp_chir.atom_id_centre", []),
                     values.get("_chem_comp_chir.volume_sign", [])))
    unsigned = [atom for atom in _SUGAR_CHIRAL_CENTRES
                if signs.get(atom) not in {"positive", "negative"}]
    if unsigned:
        raise ValueError(f"{code} must restrain the sugar chirality at {', '.join(unsigned)}")


def validate_curated_dictionary(code: str, path: Path) -> None:
    """Reject a packaged dictionary whose reviewed local topology changed."""
    try:
        ligand = CURATED_LIGANDS[code]
    except KeyError as exc:
        raise KeyError(f"No curated NASolve ligand is registered for {code}") from exc
    if ligand.reviewed_restraints:
        _validate_reviewed_restraints(ligand, path)
    bonds = _dictionary_bonds(path)
    missing = [pair for pair in ligand.required_bonds if frozenset(pair) not in bonds]
    forbidden = [pair for pair in ligand.forbidden_bonds if frozenset(pair) in bonds]
    if missing or forbidden:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join("-".join(pair) for pair in missing))
        if forbidden:
            details.append("forbidden " + ", ".join("-".join(pair) for pair in forbidden))
        raise ValueError(
            f"Curated dictionary {code} failed topology validation: {'; '.join(details)}"
        )


__all__ = [
    "CURATED_LIGAND_CODES",
    "CURATED_LIGANDS",
    "AtomSubstitution",
    "CuratedLigand",
    "RingSubstituent",
    "curated_dictionary",
    "dictionary_atom_names",
    "dictionary_ideal_bond_length",
    "ligand_data_directory",
    "validate_curated_dictionary",
]
