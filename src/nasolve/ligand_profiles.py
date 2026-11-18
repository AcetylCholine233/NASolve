"""Reviewed monomer profiles and per-site linked-phosphate modifications.

The reviewed profiles are the curated ligands with ``reviewed_restraints``
(1AP, DYI). Do not infer new chemistry for other codes. Numerical dictionary
restraints are preserved; models are never edited here.
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO

from .model_assessment import file_sha256
from .phosphate import (PhosphateError, _atoms, audit_phosphates,
                        requested_op3_sites, validate_op3_sites)
from .run_context import artifact_reference
from .backbone import requested_backbone_policy

from .curated_ligands import CURATED_LIGANDS

#: Curated ligands whose reviewed dictionary replaces ReadySet's.
AUTHORITATIVE_CODES = frozenset(
    code for code, ligand in CURATED_LIGANDS.items() if ligand.reviewed_restraints
)
MOD_ID = "NASnoOP3"
MOD_CIF = """data_mod_NASnoOP3
loop_
_chem_mod_atom.mod_id
_chem_mod_atom.function
_chem_mod_atom.atom_id
_chem_mod_atom.new_atom_id
_chem_mod_atom.new_type_symbol
_chem_mod_atom.new_type_energy
_chem_mod_atom.new_partial_charge
NASnoOP3 delete OP3 . . . .
"""


@dataclass
class _Block:
    text: str
    values: dict


def _blocks(path: Path) -> list[_Block]:
    """Split CIF1 data blocks without splitting semicolon-delimited text.

    Bio.PDB performs the token/loop parsing of each block. Do not flatten a
    multi-component file with MMCIF2Dict: repeated category names lose data.
    """
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    starts = []
    in_text = False
    for i, line in enumerate(lines):
        if line.startswith(";"):
            in_text = not in_text
        elif not in_text and re.match(r"^\s*data_\S+\s*(?:#.*)?$", line, re.I):
            starts.append(i)
    if in_text or not starts:
        raise PhosphateError(f"Invalid or unsupported dictionary CIF: {path.name}")
    result = []
    for first, last in zip(starts, [*starts[1:], len(lines)]):
        text = "".join(lines[first:last])
        try:
            values = dict(MMCIF2Dict(io.StringIO(text)))
        except (ValueError, StopIteration, IndexError) as exc:
            raise PhosphateError(f"Cannot parse dictionary {path.name}: {exc}") from exc
        result.append(_Block(text, values))
    if len({block.values["data_"].lower() for block in result}) != len(result):
        raise PhosphateError(f"Duplicate CIF data blocks in {path.name}")
    return result


def _emit(values: dict) -> str:
    output = io.StringIO()
    writer = MMCIFIO()
    writer.set_dict(values)
    writer.save(output)
    return output.getvalue()


def _component_ids(block: _Block) -> set[str]:
    return set(block.values.get("_chem_comp_atom.comp_id", []))


def normalize_ccp4_torsion_alternates(source: Path, output: Path) -> bool:
    """Materialize a Phenix-compatible copy of CCP4 alternative torsions.

    Some CCP4 monomer dictionaries encode alternative conformations as
    multiple torsion rows over the same four atoms. Phenix expects those
    alternatives folded into ``_chem_comp_tor.alt_value_angle`` instead of
    simultaneous conflicting proxies. The reviewed source is never modified.
    """
    rendered: list[str] = []
    changed = False

    atom_keys = tuple(f"_chem_comp_tor.atom_id_{index}" for index in range(1, 5))
    angle_key = "_chem_comp_tor.value_angle"
    alt_key = "_chem_comp_tor.alt_value_angle"

    for block in _blocks(source):
        values = dict(block.values)

        if alt_key in values or angle_key not in values:
            rendered.append(block.text)
            continue

        columns = [values.get(key, []) for key in atom_keys]
        angles = values.get(angle_key, [])
        count = len(angles)

        if not count or any(len(column) != count for column in columns):
            rendered.append(block.text)
            continue

        first_for_atoms: dict[tuple[str, ...], int] = {}
        keep: list[int] = []
        alternatives: dict[int, list[str]] = {}

        for index in range(count):
            atom_set = tuple(sorted(column[index] for column in columns))
            first = first_for_atoms.get(atom_set)
            if first is None:
                first_for_atoms[atom_set] = index
                keep.append(index)
                alternatives[index] = []
                continue

            changed = True
            angle = angles[index]
            if angle != angles[first] and angle not in alternatives[first]:
                alternatives[first].append(angle)

        if len(keep) == count:
            rendered.append(block.text)
            continue

        for key, column in list(values.items()):
            if key.startswith("_chem_comp_tor.") and isinstance(column, list):
                values[key] = [column[index] for index in keep]

        values[alt_key] = [
            ",".join(alternatives[index]) if alternatives[index] else "."
            for index in keep
        ]
        rendered.append(_emit(values))

    with output.open("x", encoding="utf-8") as handle:
        handle.write(
            "# Runtime Phenix compatibility adaptation: CCP4 duplicate torsion "
            "targets folded into _chem_comp_tor.alt_value_angle. "
            "Reviewed source unchanged.\n"
        )
        handle.write("\n".join(rendered))

    return changed


def combine_dictionary_inputs(paths: Sequence[Path], output: Path) -> None:
    """Merge component lists; never concatenate duplicate data_comp_list blocks."""
    blocks = []
    names = set()
    list_rows: list[dict[str, str]] = []
    keys = []
    seen_ids = set()
    for path in paths:
        for block in _blocks(path):
            values = block.values
            if "_chem_comp.id" in values and not _component_ids(block):
                columns = {k: v for k, v in values.items() if k.startswith("_chem_comp.")}
                n = len(values["_chem_comp.id"])
                if any(len(v) != n for v in columns.values()):
                    raise PhosphateError(f"Inconsistent component list in {path.name}")
                for key in columns:
                    if key not in keys:
                        keys.append(key)
                for i in range(n):
                    row = {k: v[i] for k, v in columns.items()}
                    if row["_chem_comp.id"] in seen_ids:
                        raise PhosphateError("Duplicate component identity in supplied dictionaries")
                    seen_ids.add(row["_chem_comp.id"])
                    list_rows.append(row)
            else:
                name = values["data_"].lower()
                if name in names:
                    raise PhosphateError(f"Conflicting dictionary data block: {name}")
                names.add(name)
                blocks.append(block.text)
    prefix = ""
    if list_rows:
        prefix = _emit({"data_": "comp_list", **{
            k: [row.get(k, "?") for row in list_rows] for k in keys
        }})
    with output.open("x", encoding="utf-8") as handle:
        handle.write(prefix + "\n" + "\n".join(blocks))


def _exclude_components(source: Path, codes: set[str], output: Path) -> tuple[Path | None, set[str]]:
    """Keep ReadySet additions except authoritative component definitions."""
    retained = []
    present = set()
    for block in _blocks(source):
        values = block.values
        ids = _component_ids(block)
        if ids & codes:
            if ids - codes:
                raise PhosphateError("Mixed-component atom block cannot be overridden safely")
            continue
        if not ids and "_chem_comp.id" in values:
            keep = [i for i, code in enumerate(values["_chem_comp.id"]) if code not in codes]
            if not keep:
                continue
            if len(keep) != len(values["_chem_comp.id"]):
                if any(k != "data_" and not k.startswith("_chem_comp.") for k in values):
                    raise PhosphateError("Mixed-category component list requires review")
                values = {k: ([v[i] for i in keep] if k != "data_" else v)
                          for k, v in values.items()}
                retained.append(_emit(values))
                continue
        retained.append(block.text)
        present.update(ids)
    if not retained:
        return None, present
    with output.open("x", encoding="utf-8") as handle:
        handle.write("# ReadySet output with reviewed component overrides excluded.\n" + "\n".join(retained))
    return output, present


def _profile_inventory(model: Path) -> dict[str, str]:
    residues = {}
    for atom in _atoms(model.read_text(encoding="utf-8").splitlines(keepends=True)):
        code = atom.line[17:20].strip()
        if code in AUTHORITATIVE_CODES:
            residues[atom.site] = code
    return residues


def write_linked_profile(
    model: Path, directory: Path, *, allow_op3_sites: tuple[str, ...] = (),
    passthrough_sites: tuple[str, ...] = (),
) -> tuple[dict[str, object], tuple[Path, ...]]:
    inventory = _profile_inventory(model)
    audit = audit_phosphates(
        model, allow_op3_sites=allow_op3_sites, inspect_sites=tuple(inventory),
        passthrough_sites=passthrough_sites,
    )
    modifications = []
    for item in audit["checked"]:
        if item["site"] in inventory and item["incoming_site"] is not None:
            chain, resid = item["site"].split(":")
            # Supported selectors are explicit PDB identities, not arbitrary PHIL.
            validate_op3_sites([item["site"]])
            selection = f"chain '{'' if chain == '_' else chain}' and resid {resid} and resname {inventory[item['site']]}"
            modifications.append({"site": item["site"], "residue": inventory[item["site"]],
                                  "incoming_site": item["incoming_site"], "outgoing_site": item["outgoing_site"], "data_mod": MOD_ID,
                                  "selection": selection})
    paths = ()
    if modifications:
        cif = directory / "linked_phosphate.cif"
        phil = directory / "linked_phosphate.phil"
        with cif.open("x", encoding="utf-8") as handle:
            handle.write(MOD_CIF)
        with phil.open("x", encoding="utf-8") as handle:
            handle.write("# Generated from verified connectivity; do not apply by residue code alone.\n")
            for item in modifications:
                handle.write("refinement.pdb_interpretation.apply_cif_modification {\n"
                             f"  data_mod = {MOD_ID}\n  residue_selection = {item['selection']}\n}}\n")
        paths = (cif, phil)
    return {"schema_version": 1, "inventory": inventory, "modifications": modifications}, paths


def validate_model_phosphate_policy(model: Path, report: Mapping[str, object]) -> dict[str, object]:
    """Validate both default OP3 policy and the identities frozen in the PHIL."""
    allowed = requested_op3_sites(report)
    backbone = requested_backbone_policy(report)
    passthrough = tuple(backbone["experimental_passthrough_sites"])
    postmr = report.get("postmr")
    profile = postmr.get("linked_phosphate_profile") if isinstance(postmr, Mapping) else None
    modifications = []
    if profile is not None:
        if not isinstance(profile, Mapping) or profile.get("schema_version") != 1:
            raise PhosphateError("Malformed linked-phosphate profile")
        if profile.get("inventory") != _profile_inventory(model):
            raise PhosphateError("Reviewed nucleotide identities changed; prepare new dictionary profiles")
        modifications = profile.get("modifications")
        if not isinstance(modifications, list) or any(
            not isinstance(item, Mapping) or not isinstance(item.get("site"), str)
            or not isinstance(item.get("incoming_site"), str) or item.get("data_mod") != MOD_ID
            for item in modifications
        ):
            raise PhosphateError("Malformed linked-phosphate modification inventory")
    result = audit_phosphates(
        model, allow_op3_sites=allowed,
        check_sites=tuple(item["site"] for item in modifications),
        passthrough_sites=passthrough,
    )
    actual = {item["site"]: (item["incoming_site"], item["outgoing_site"]) for item in result["checked"]}
    if any(actual.get(item["site"]) != (item["incoming_site"], item.get("outgoing_site"))
           for item in modifications):
        raise PhosphateError("Phosphate backbone link differs from the frozen site modification")
    return result


def effective_restraints(
    paths: Sequence[Path], generated: Path | None, directory: Path,
) -> tuple[list[Path], list[Path]]:
    """Resolve CIF precedence once, before checkpoint creation (not at reuse)."""
    cifs = [p for p in paths if p.suffix.lower() == ".cif" and p.name != "linked_phosphate.cif"]
    codes = {p.stem for p in cifs} & AUTHORITATIVE_CODES
    other = [p for p in paths if p not in cifs]
    dictionaries = []
    generated_codes = set()
    if generated is not None:
        filtered, generated_codes = _exclude_components(
            generated, codes, directory / "readyset_supplement.cif")
        if filtered is not None:
            dictionaries.append(filtered)
    for cif in cifs:
        if cif.stem in codes or generated is None or cif.stem not in generated_codes:
            dictionaries.append(cif)
    return [*other, *dictionaries], dictionaries


def frozen_reference(path: Path, run: Path) -> dict[str, object]:
    return {**artifact_reference(path, run), "sha256": file_sha256(path), "size": path.stat().st_size}


def is_modification_only(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    return "_chem_mod_atom." in text and "_chem_comp_atom." not in text
