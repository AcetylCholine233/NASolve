#!/usr/bin/env python3
"""Build NASolve's DYI restraint dictionary (O4-triazolyl-dU) from AceDRG entries.

DYI is the nucleotide suppliers sell as O4-triazolyl-dU (a "convertible"
nucleoside): 4-(1,2,4-triazol-1-yl)-1-(2-deoxy-5-O-phosphono-beta-D-erythro-
pentofuranosyl)pyrimidin-2(1H)-one. Its ring has dC's bonding (N3=C4, a
nitrogen on C4), so it is built as dC with the N4 amino group replaced by a
1,2,4-triazol-1-yl ring whose N1 is DYI's N4. DYI is not a wwPDB CCD
component (checked against the CCD in Phenix 2.2.1-6174, entries to
2026-07-10); DYI is a laboratory label.

Sources, both from the CCP4 Monomer Library (AceDRG dictionaries):

- ``d/DC.cif``: every atom shared with dC (sugar, phosphate, pyrimidinone
  ring, H5/H6) with its bonds, angles, torsions, chiral centres and base
  plane, copied unchanged. Only N4 changes: energy type NR5 (as TZF's N1) and
  its two hydrogens are dropped.
- ``t/TZF.cif`` (4-(1,2,4-triazol-1-yl)phenol): the triazole ring, its two
  hydrogens, the angles at N1 including those to the ring carbon, the ring
  torsions, the triazole plane and the ring-to-ring torsion (coplanar,
  period 2).

One value is chosen here: C4-N4 = 1.411 +/- 0.015 A. AceDRG gives
1.424 A for a triazol-1-yl on a carbon between two carbons (TZF, ZRS; 1.429 A
on quinolin-4-yl, TQU) and 1.399 A on a carbon between two nitrogens
(triazinyl, 1S1). DYI's C4 sits between N3 and C5, so the mean of the two
classes is used with an e.s.d. covering both.

Coordinates: DC's AceDRG coordinates for every shared atom; N4 moved along
dC's C4->N4 line to 1.411 A; the triazole atoms and hydrogens fitted to the
restraints by least squares, coplanar with the base and with N2T on the C5
side (C5-C4-N4-N2T = 0). The other coplanar orientation puts H5T against H5
and the N2T and N3 lone pairs together.

Usage::

    python scripts/build_dyi_dictionary.py DC.cif TZF.cif src/nasolve/data/ligands/DYI.cif
"""

from __future__ import annotations

import argparse
import hashlib
import io
import math
import re
import sys
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.optimize import least_squares

CODE = "DYI"
NAME = "O4-TRIAZOLYL-2'-DEOXYURIDINE-5'-MONOPHOSPHATE"
C4_N4 = (1.411, 0.015)
#: TZF atom -> DYI atom. CAK (the phenol ring carbon) is DYI's C4.
TZF_MAP = {"NAL": "N4", "NAI": "N2T", "CAF": "C3T", "NAH": "N4T", "CAG": "C5T",
           "H6": "H5T", "H7": "H3T", "CAK": "C4"}
TRIAZOLE = ("N2T", "C3T", "N4T", "C5T")
MOVING = ("N4", *TRIAZOLE, "H3T", "H5T")
DROPPED = {"H41", "H42"}
#: The one TZF torsion that also names a phenol atom: CAE CAK NAL CAG (ortho
#: carbon, ring carbon, N1, C5). DYI's ortho carbon is C5.
INTER_RING = ("C5", "C4", "N4", "C5T")
SMILES = "n3cncn3C1=NC(=O)N(C=C1)[C@H]2C[C@H](O)[C@@H](CO[P](O)(O)=O)O2"


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def component_block(path: Path, code: str) -> dict:
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"^data_comp_{re.escape(code)}\s*$", text, re.M)
    if match is None:
        raise SystemExit(f"{path} has no data_comp_{code} block")
    following = re.search(r"^data_", text[match.end():], re.M)
    block = text[match.start(): match.end() + following.start() if following else None]
    return dict(MMCIF2Dict(io.StringIO(block)))


def rows(block: dict, category: str, fields: tuple[str, ...]) -> list[dict]:
    columns = [block[f"{category}.{field}"] for field in fields]
    return [dict(zip(fields, values)) for values in zip(*columns)]


ATOM = ("comp_id", "atom_id", "type_symbol", "type_energy", "charge", "x", "y", "z")
BOND = ("comp_id", "atom_id_1", "atom_id_2", "type", "aromatic", "value_dist_nucleus",
        "value_dist_nucleus_esd", "value_dist", "value_dist_esd")
ANGLE = ("comp_id", "atom_id_1", "atom_id_2", "atom_id_3", "value_angle", "value_angle_esd")
TOR = ("comp_id", "id", "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
       "value_angle", "value_angle_esd", "period")
CHIR = ("comp_id", "id", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3", "volume_sign")
PLANE = ("comp_id", "plane_id", "atom_id", "dist_esd")


def assemble(dc: dict, tzf: dict) -> dict[str, list[dict]]:
    names = lambda row, keys: [row[key] for key in keys]  # noqa: E731
    out: dict[str, list[dict]] = {}

    atoms = []
    for row in rows(dc, "_chem_comp_atom", ATOM):
        if row["atom_id"] in DROPPED:
            continue
        row = {**row, "comp_id": CODE}
        if row["atom_id"] == "N4":
            row["type_energy"] = "NR5"
        atoms.append(row)
        if row["atom_id"] == "C6":   # triazole heavy atoms after the base
            atoms.extend({"comp_id": CODE, "atom_id": name, "type_symbol": name[0],
                          "type_energy": None, "charge": "0", "x": None, "y": None, "z": None}
                         for name in TRIAZOLE)
    atoms.extend({"comp_id": CODE, "atom_id": name, "type_symbol": "H", "type_energy": "H",
                  "charge": "0", "x": None, "y": None, "z": None} for name in ("H3T", "H5T"))
    tzf_types = {TZF_MAP[row["atom_id"]]: row["type_energy"]
                 for row in rows(tzf, "_chem_comp_atom", ATOM) if row["atom_id"] in TZF_MAP}
    for row in atoms:
        if row["atom_id"] in TRIAZOLE:
            row["type_energy"] = tzf_types[row["atom_id"]]
    out["atom"] = atoms

    bonds = []
    for row in rows(dc, "_chem_comp_bond", BOND):
        if DROPPED & set(names(row, ("atom_id_1", "atom_id_2"))):
            continue
        row = {**row, "comp_id": CODE}
        if {row["atom_id_1"], row["atom_id_2"]} == {"C4", "N4"}:
            value, esd = C4_N4
            row.update(value_dist_nucleus=f"{value:.3f}", value_dist_nucleus_esd=f"{esd:.4f}",
                       value_dist=f"{value:.3f}", value_dist_esd=f"{esd:.4f}")
        bonds.append(row)
    for row in rows(tzf, "_chem_comp_bond", BOND):
        pair = names(row, ("atom_id_1", "atom_id_2"))
        if all(atom in TZF_MAP for atom in pair) and "CAK" not in pair:
            bonds.append({**row, "comp_id": CODE, "atom_id_1": TZF_MAP[pair[0]],
                          "atom_id_2": TZF_MAP[pair[1]]})
    out["bond"] = bonds

    angles = [{**row, "comp_id": CODE} for row in rows(dc, "_chem_comp_angle", ANGLE)
              if not DROPPED & set(names(row, ("atom_id_1", "atom_id_2", "atom_id_3")))]
    for row in rows(tzf, "_chem_comp_angle", ANGLE):
        triple = names(row, ("atom_id_1", "atom_id_2", "atom_id_3"))
        if all(atom in TZF_MAP for atom in triple) and triple[1] != "CAK":
            angles.append({**row, "comp_id": CODE,
                           **{f"atom_id_{i + 1}": TZF_MAP[atom] for i, atom in enumerate(triple)}})
    out["angle"] = angles

    torsions = [{**row, "comp_id": CODE} for row in rows(dc, "_chem_comp_tor", TOR)
                if not DROPPED & set(names(row, ("atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4")))]
    for row in rows(tzf, "_chem_comp_tor", TOR):
        quad = names(row, ("atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4"))
        if quad == ["CAE", "CAK", "NAL", "CAG"]:
            mapped = list(INTER_RING)
        elif all(atom in TZF_MAP for atom in quad):
            mapped = [TZF_MAP[atom] for atom in quad]
        else:
            continue
        torsions.append({**row, "comp_id": CODE, "id": f"tz_{row['id']}",
                         **{f"atom_id_{i + 1}": atom for i, atom in enumerate(mapped)}})
    out["tor"] = torsions

    out["chir"] = [{**row, "comp_id": CODE} for row in rows(dc, "_chem_comp_chir", CHIR)]

    planes = []
    dc_planes = rows(dc, "_chem_comp_plane_atom", PLANE)
    for plane_id in dict.fromkeys(row["plane_id"] for row in dc_planes):
        members = [row for row in dc_planes if row["plane_id"] == plane_id]
        if not DROPPED & {row["atom_id"] for row in members}:   # drop dC's NH2 plane
            planes.extend({**row, "comp_id": CODE} for row in members)
    tzf_planes = rows(tzf, "_chem_comp_plane_atom", PLANE)
    triazole_plane = [row for row in tzf_planes if row["atom_id"] == "NAI"][0]["plane_id"]
    planes.extend({**row, "comp_id": CODE, "plane_id": "plan-2", "atom_id": TZF_MAP[row["atom_id"]]}
                  for row in tzf_planes if row["plane_id"] == triazole_plane)
    out["plane"] = planes
    return out


# ---------------------------------------------------------------- geometry

def dihedral(a, b, c, d) -> float:
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def angle(a, b, c) -> float:
    u, v = a - b, c - b
    return math.degrees(math.acos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))


def plane_deviations(points: np.ndarray) -> np.ndarray:
    centred = points - points.mean(axis=0)
    normal = np.linalg.svd(centred)[2][-1]
    return centred @ normal


def wrapped(value: float) -> float:
    return (value + 180.0) % 360.0 - 180.0


def place_triazole(xyz: dict[str, np.ndarray], restraints: dict[str, list[dict]]) -> None:
    """Start the triazole in the base plane, then fit it to the restraints."""
    c4 = xyz["C4"]
    u = (xyz["N4"] - c4) / np.linalg.norm(xyz["N4"] - c4)
    base = np.array([xyz[name] for name in ("N1", "C2", "N3", "C4", "C5", "C6")])
    normal = np.linalg.svd(base - base.mean(axis=0))[2][-1]
    v = np.cross(normal, u)
    if np.dot(xyz["C5"] - c4, v) < 0:   # v points to the C5 side
        v = -v
    n4 = c4 + C4_N4[0] * u
    xyz["N4"] = n4

    def towards(origin, back, degrees, side, length):
        back = (back - origin) / np.linalg.norm(back - origin)
        perp = side - np.dot(side, back) * back
        perp /= np.linalg.norm(perp)
        theta = math.radians(degrees)
        return origin + length * (math.cos(theta) * back + math.sin(theta) * perp)

    xyz["N2T"] = towards(n4, c4, 121.5, v, 1.362)      # N2T on the C5 side
    xyz["C5T"] = towards(n4, c4, 129.2, -v, 1.335)
    xyz["C3T"] = towards(xyz["N2T"], n4, 101.8, -v, 1.313)
    xyz["N4T"] = towards(xyz["C5T"], n4, 110.6, v, 1.314)
    xyz["H3T"] = xyz["C3T"] + 0.936 * (xyz["C3T"] - (xyz["N2T"] + xyz["N4T"]) / 2) / np.linalg.norm(
        xyz["C3T"] - (xyz["N2T"] + xyz["N4T"]) / 2)
    xyz["H5T"] = xyz["C5T"] + 0.942 * (xyz["C5T"] - (n4 + xyz["N4T"]) / 2) / np.linalg.norm(
        xyz["C5T"] - (n4 + xyz["N4T"]) / 2)

    moving = set(MOVING)
    bonds = [(r["atom_id_1"], r["atom_id_2"], float(r["value_dist"]), float(r["value_dist_esd"]))
             for r in restraints["bond"] if moving & {r["atom_id_1"], r["atom_id_2"]}]
    angles = [(r["atom_id_1"], r["atom_id_2"], r["atom_id_3"], float(r["value_angle"]),
               float(r["value_angle_esd"]))
              for r in restraints["angle"] if moving & {r["atom_id_1"], r["atom_id_2"], r["atom_id_3"]}]
    torsions = [(tuple(r[f"atom_id_{i}"] for i in range(1, 5)), float(r["value_angle"]),
                 float(r["value_angle_esd"]), int(r["period"]))
                for r in restraints["tor"]
                if moving & {r[f"atom_id_{i}"] for i in range(1, 5)}]
    planes = {}
    for r in restraints["plane"]:
        planes.setdefault(r["plane_id"], []).append((r["atom_id"], float(r["dist_esd"])))
    planes = {pid: members for pid, members in planes.items() if moving & {m for m, _ in members}}

    def coordinates(vector):
        current = dict(xyz)
        for index, name in enumerate(MOVING):
            current[name] = vector[3 * index: 3 * index + 3]
        return current

    def residuals(vector):
        current = coordinates(vector)
        out = [(np.linalg.norm(current[a] - current[b]) - d) / s for a, b, d, s in bonds]
        out += [(angle(current[a], current[b], current[c]) - t) / s for a, b, c, t, s in angles]
        out += [wrapped(period * (dihedral(*(current[n] for n in quad)) - t)) / period / s
                for quad, t, s, period in torsions]
        for members in planes.values():
            points = np.array([current[name] for name, _ in members])
            out += list(plane_deviations(points) / np.array([s for _, s in members]))
        return np.array(out)

    start = np.concatenate([xyz[name] for name in MOVING])
    fit = least_squares(residuals, start, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    xyz.update({name: value for name, value in coordinates(fit.x).items() if name in MOVING})


# ---------------------------------------------------------------- output

def quoted(name: str) -> str:
    return f'"{name}"' if "'" in name else name


def table(category: str, fields: tuple[str, ...], data: list[dict]) -> str:
    cells = [[quoted(str(row[field])) for field in fields] for row in data]
    widths = [max(len(row[i]) for row in cells) for i in range(len(fields))]
    lines = ["loop_", *(f"{category}.{field}" for field in fields)]
    lines += [" ".join(cell.ljust(width) for cell, width in zip(row, widths)).rstrip() for row in cells]
    return "\n".join(lines)


def write(path: Path, restraints: dict[str, list[dict]], header: list[str]) -> None:
    heavy = sum(1 for row in restraints["atom"] if row["type_symbol"] != "H")
    parts = [
        "\n".join(f"# {line}" if line else "#" for line in header),
        "data_comp_list",
        table("_chem_comp", ("id", "three_letter_code", "name", "group", "number_atoms_all",
                             "number_atoms_nh", "desc_level"),
              [{"id": CODE, "three_letter_code": CODE, "name": f'"{NAME}"', "group": "DNA",
                "number_atoms_all": len(restraints["atom"]), "number_atoms_nh": heavy,
                "desc_level": "."}]).replace(f'""{NAME}""', f'"{NAME}"'),
        "",
        f"data_comp_{CODE}",
        table("_chem_comp_atom", ATOM, restraints["atom"]),
        table("_chem_comp_bond", BOND, restraints["bond"]),
        table("_chem_comp_angle", ANGLE, restraints["angle"]),
        table("_chem_comp_tor", TOR, restraints["tor"]),
        table("_chem_comp_chir", CHIR, restraints["chir"]),
        table("_chem_comp_plane_atom", PLANE, restraints["plane"]),
        table("_pdbx_chem_comp_descriptor", ("comp_id", "type", "program", "program_version",
                                            "descriptor"),
              [{"comp_id": CODE, "type": "SMILES", "program": "NASolve",
                "program_version": "1", "descriptor": SMILES}]),
    ]
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dc", type=Path, help="CCP4 Monomer Library d/DC.cif")
    parser.add_argument("tzf", type=Path, help="CCP4 Monomer Library t/TZF.cif")
    parser.add_argument("output", type=Path)
    parser.add_argument("--library-commit", default="unknown",
                        help="Monomer Library commit the two files come from")
    args = parser.parse_args(argv)
    dc, tzf = component_block(args.dc, "DC"), component_block(args.tzf, "TZF")
    restraints = assemble(dc, tzf)
    xyz = {row["atom_id"]: np.array([float(row["x"]), float(row["y"]), float(row["z"])])
           for row in restraints["atom"] if row["x"] is not None}
    place_triazole(xyz, restraints)
    for row in restraints["atom"]:
        row["x"], row["y"], row["z"] = (f"{value:.3f}" for value in xyz[row["atom_id"]])
    header = [
        f"NASolve reviewed restraint dictionary {CODE}: O4-triazolyl-dU (supplier name),",
        "4-(1,2,4-triazol-1-yl)-1-(2-deoxy-5-O-phosphono-beta-D-erythro-pentofuranosyl)pyrimidin-2(1H)-one.",
        "Not a wwPDB CCD component; DYI is a laboratory label. Parent DC (N4 = triazole N1).",
        f"Built by scripts/build_dyi_dictionary.py from the CCP4 Monomer Library (AceDRG), commit {args.library_commit}:",
        f"  d/DC.cif  blob {git_blob_sha1(args.dc)}: all atoms shared with dC, unchanged",
        f"  t/TZF.cif blob {git_blob_sha1(args.tzf)}: triazol-1-yl ring, angles at N1, ring torsions and plane",
        f"Chosen here: C4-N4 {C4_N4[0]:.3f} +/- {C4_N4[1]:.3f} A (AceDRG 1.424 on C(C,C), 1.399 on C(N,N));",
        "N4 energy type NR5. Coordinates: DC's for shared atoms; triazole fitted, coplanar, N2T on the C5 side.",
    ]
    write(args.output, restraints, header)
    print(f"wrote {args.output}: {len(restraints['atom'])} atoms, {len(restraints['bond'])} bonds, "
          f"{len(restraints['angle'])} angles, {len(restraints['tor'])} torsions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
