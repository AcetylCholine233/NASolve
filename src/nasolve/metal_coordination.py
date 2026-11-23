"""Metal coordination restraints for AutoRefine (``--coordination``).

``nasolve autorefine RUN --coordination "B:101=B:4/N2T,A:12/O4"`` restrains a
heavy-atom site to its binding partners (donors) in phenix.refine's own form,
``refinement.geometry_restraints.edits``:

* one ``bond`` per donor, its ideal length from :data:`COORDINATION_LENGTHS`
  (by metal and donor element), sigma :data:`BOND_SIGMA`;
* for two donors (linear coordination), one ``angle`` donor-donor-metal with
  the first donor at the vertex, ideal 0 degrees, sigma 3. It is not a
  physical angle: at 0 degrees the metal lies on the line through the two
  donors, which keeps the donor-metal-donor geometry linear without the
  180-degree angle at the metal;
* for three or more donors, bonds only (no preset angle).

Values chosen by the lab for now (user, 2026-11-23). Phenix adds custom bonds
to its bonding list, so the metal and its donors are no longer treated as
clashing. Phenix refuses ``action = add`` for a bond it already has (metals
are not bonded automatically with ``link_metals = Auto``, AutoRefine's
setting).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

#: Ideal metal-donor bond lengths (A) by (metal element, donor element).
#: Lab values for now (user, 2026-11-23).
COORDINATION_LENGTHS: dict[tuple[str, str], float] = {
    ("AG", "O"): 2.5,
    ("AG", "N"): 2.3,
    ("HG", "O"): 2.1,
    ("HG", "N"): 2.1,
}
BOND_SIGMA = 0.2
#: Linear coordination: angle donor-donor-metal (vertex at the first donor).
LINEAR_ANGLE_IDEAL = 0.0
LINEAR_ANGLE_SIGMA = 3.0
#: A donor this much farther from the metal than its ideal length is reported.
FAR_DONOR_MARGIN = 1.5

_SITE = re.compile(r"^\s*([A-Za-z0-9]{1,4}):(-?\d+[A-Za-z]?)(?:/([^,=\s/]+))?\s*$")


class CoordinationError(ValueError):
    """A --coordination specification cannot be used."""


@dataclass(frozen=True)
class AtomRef:
    chain: str
    resid: str
    name: str | None    # None: the metal residue's only atom

    @property
    def label(self) -> str:
        return f"{self.chain}:{self.resid}" + (f"/{self.name}" if self.name else "")


@dataclass(frozen=True)
class CoordinationSpec:
    text: str
    metal: AtomRef
    donors: tuple[AtomRef, ...]


@dataclass(frozen=True)
class ModelAtom:
    chain: str
    resid: str
    resname: str
    name: str
    element: str
    xyz: tuple[float, float, float]

    @property
    def label(self) -> str:
        return f"{self.resname} {self.chain}:{self.resid} {self.name}"

    @property
    def selection(self) -> str:
        return f"chain {self.chain} and resid {self.resid} and name {self.name}"


@dataclass(frozen=True)
class CoordinationBond:
    metal: ModelAtom
    donor: ModelAtom
    ideal: float
    sigma: float
    current: float


@dataclass(frozen=True)
class CoordinationAngle:
    atoms: tuple[ModelAtom, ModelAtom, ModelAtom]   # donor 2, donor 1 (vertex), metal
    ideal: float
    sigma: float
    current: float


@dataclass(frozen=True)
class CoordinationSite:
    spec: CoordinationSpec
    metal: ModelAtom
    bonds: tuple[CoordinationBond, ...]
    angle: CoordinationAngle | None
    warnings: tuple[str, ...]

    @property
    def geometry(self) -> str:
        count = len(self.bonds)
        if count == 2:
            return "linear"
        return "one donor" if count == 1 else f"{count} donors (no angle restraint)"


def _table_text() -> str:
    return ", ".join(f"{metal.title()}-{donor} {length:g} A"
                     for (metal, donor), length in COORDINATION_LENGTHS.items())


def _site(text: str, what: str) -> AtomRef:
    match = _SITE.match(text)
    if not match:
        raise CoordinationError(
            f"{what} {text.strip()!r} is not CHAIN:RESID or CHAIN:RESID/ATOM (e.g. B:4/N2T)"
        )
    chain, resid, name = match.groups()
    return AtomRef(chain, resid, name)


def parse_coordination(text: str) -> CoordinationSpec:
    """``METAL=DONOR,DONOR``: e.g. ``B:101=B:4/N2T,A:12/O4``. The metal's atom
    name may be left out when its residue has one atom; donors need theirs."""
    if text.count("=") != 1:
        raise CoordinationError(
            f"--coordination {text!r}: give METAL=DONOR[,DONOR...], e.g. B:101=B:4/N2T,A:12/O4"
        )
    metal_text, donor_text = text.split("=")
    metal = _site(metal_text, "Metal")
    donors = tuple(_site(part, "Donor") for part in donor_text.split(",") if part.strip())
    if not donors:
        raise CoordinationError(f"--coordination {text!r}: no donor atoms given")
    for donor in donors:
        if donor.name is None:
            raise CoordinationError(
                f"--coordination {text!r}: donor {donor.label} needs its atom name (e.g. {donor.label}/N3)"
            )
    if len(set(donors)) != len(donors):
        raise CoordinationError(f"--coordination {text!r}: a donor is given twice")
    return CoordinationSpec(text=text.strip(), metal=metal, donors=donors)


def _element(line: str) -> str:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    name = line[12:16]
    return (name[:2] if name[0] != " " else name[1:2]).strip().upper()


def read_model_atoms(model: Path) -> dict[tuple[str, str], list[tuple[str, ModelAtom]]]:
    """Atoms of the first model by (chain, residue number): (altloc, atom)."""
    residues: dict[tuple[str, str], list[tuple[str, ModelAtom]]] = {}
    for line in model.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        chain, resid = line[21].strip(), line[22:27].strip()
        atom = ModelAtom(
            chain=chain, resid=resid, resname=line[17:20].strip(), name=line[12:16].strip(),
            element=_element(line),
            xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
        )
        residues.setdefault((chain, resid), []).append((line[16].strip(), atom))
    return residues


def _find(residues, ref: AtomRef, what: str) -> ModelAtom:
    atoms = residues.get((ref.chain, ref.resid))
    if not atoms:
        raise CoordinationError(f"{what} {ref.label}: no residue {ref.chain}:{ref.resid} in the model")
    if ref.name is None:
        names = sorted({atom.name for _, atom in atoms})
        if len(names) != 1:
            raise CoordinationError(
                f"{what} {ref.label}: the residue has several atoms ({', '.join(names)}); "
                f"name one, e.g. {ref.label}/{names[0]}"
            )
        name = names[0]
    else:
        name = ref.name
    if "'" in name or '"' in name:
        raise CoordinationError(f"{what} {ref.label}: atom names with quotes are not supported yet")
    found = [(altloc, atom) for altloc, atom in atoms if atom.name == name]
    if not found:
        raise CoordinationError(
            f"{what} {ref.label}: no atom {name} in {atoms[0][1].resname} {ref.chain}:{ref.resid} "
            f"(atoms: {', '.join(sorted({atom.name for _, atom in atoms}))})"
        )
    if len(found) > 1:
        raise CoordinationError(
            f"{what} {ref.label}: has alternate conformations; not supported yet"
        )
    return found[0][1]


def _angle(vertex, first, second) -> float:
    u = [a - b for a, b in zip(first, vertex)]
    v = [a - b for a, b in zip(second, vertex)]
    cosine = sum(a * b for a, b in zip(u, v)) / (math.hypot(*u) * math.hypot(*v))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def plan_coordination(model: Path, specs: Sequence[CoordinationSpec]) -> list[CoordinationSite]:
    """Resolve each specification in ``model`` and build its restraints."""
    residues = read_model_atoms(model)
    sites: list[CoordinationSite] = []
    seen: set[tuple[str, str, str]] = set()
    for spec in specs:
        metal = _find(residues, spec.metal, "Metal")
        key = (metal.chain, metal.resid, metal.name)
        if key in seen:
            raise CoordinationError(f"--coordination: {metal.label} is given twice")
        seen.add(key)
        bonds = []
        warnings = []
        for ref in spec.donors:
            donor = _find(residues, ref, "Donor")
            if (donor.chain, donor.resid, donor.name) == key:
                raise CoordinationError(f"--coordination {spec.text!r}: the metal is also a donor")
            ideal = COORDINATION_LENGTHS.get((metal.element, donor.element))
            if ideal is None:
                raise CoordinationError(
                    f"--coordination {spec.text!r}: no {metal.element.title()}-{donor.element} "
                    f"length in the table ({_table_text()})"
                )
            current = math.dist(metal.xyz, donor.xyz)
            if current > ideal + FAR_DONOR_MARGIN:
                warnings.append(
                    f"{donor.label} is {current:.2f} A from {metal.label} (ideal {ideal:g} A)"
                )
            bonds.append(CoordinationBond(metal, donor, ideal, BOND_SIGMA, round(current, 3)))
        angle = None
        if len(bonds) == 2:
            first, second = bonds[0].donor, bonds[1].donor
            angle = CoordinationAngle(
                atoms=(second, first, metal), ideal=LINEAR_ANGLE_IDEAL, sigma=LINEAR_ANGLE_SIGMA,
                current=round(_angle(first.xyz, second.xyz, metal.xyz), 2),
            )
        sites.append(CoordinationSite(spec, metal, tuple(bonds), angle, tuple(warnings)))
    return sites


def coordination_phil(sites: Sequence[CoordinationSite]) -> str:
    """The restraints as a phenix.refine parameter file."""
    lines = [
        "# NASolve metal coordination restraints (nasolve autorefine --coordination)",
        "# Bonds: lab table, sigma %g A. Linear (two donors): angle donor-donor-metal "
        "%g deg, sigma %g (not a physical angle; keeps the metal on the donor-donor line)."
        % (BOND_SIGMA, LINEAR_ANGLE_IDEAL, LINEAR_ANGLE_SIGMA),
        "refinement {",
        "  geometry_restraints.edits {",
    ]
    for site in sites:
        lines.append(f"    # {site.metal.label}: {site.spec.text} ({site.geometry})")
        for bond in site.bonds:
            lines += [
                "    bond {",
                "      action = add",
                f"      atom_selection_1 = {bond.metal.selection}",
                f"      atom_selection_2 = {bond.donor.selection}",
                f"      distance_ideal = {bond.ideal:g}",
                f"      sigma = {bond.sigma:g}",
                "    }",
            ]
        if site.angle is not None:
            first, vertex, last = site.angle.atoms
            lines += [
                "    angle {",
                "      action = add",
                f"      atom_selection_1 = {first.selection}",
                f"      atom_selection_2 = {vertex.selection}",
                f"      atom_selection_3 = {last.selection}",
                f"      angle_ideal = {site.angle.ideal:g}",
                f"      sigma = {site.angle.sigma:g}",
                "    }",
            ]
    lines += ["  }", "}", ""]
    return "\n".join(lines)


def coordination_record(site: CoordinationSite) -> dict[str, object]:
    """What a round and its checkpoint record for one site."""
    return {
        "spec": site.spec.text,
        "metal": site.metal.label,
        "metal_key": [site.metal.chain, site.metal.resid, site.metal.name],
        "geometry": site.geometry,
        "bonds": [
            {"donor": bond.donor.label, "elements": f"{bond.metal.element}-{bond.donor.element}",
             "ideal": bond.ideal, "sigma": bond.sigma, "current": bond.current}
            for bond in site.bonds
        ],
        "angle": None if site.angle is None else {
            "atoms": [atom.label for atom in site.angle.atoms],
            "ideal": site.angle.ideal, "sigma": site.angle.sigma, "current": site.angle.current,
        },
        "warnings": list(site.warnings),
        "summary": describe_site(site),
    }


def describe_site(site: CoordinationSite) -> str:
    """One line for the terminal: bonds with current and ideal lengths."""
    bonds = ", ".join(f"{bond.donor.chain}:{bond.donor.resid} {bond.donor.name} "
                      f"{bond.ideal:g} A (now {bond.current:.2f})" for bond in site.bonds)
    text = f"{site.metal.label}: {bonds}"
    if site.angle is not None:
        first, vertex, metal = site.angle.atoms
        text += (f"; linear: {first.name}-{vertex.name}-{metal.name} "
                 f"{site.angle.ideal:g} +/- {site.angle.sigma:g} deg (now {site.angle.current:.1f})")
    return text


__all__ = [
    "BOND_SIGMA", "COORDINATION_LENGTHS", "LINEAR_ANGLE_IDEAL", "LINEAR_ANGLE_SIGMA",
    "CoordinationError", "CoordinationSite", "CoordinationSpec", "coordination_phil",
    "coordination_record", "describe_site", "parse_coordination", "plan_coordination",
]
