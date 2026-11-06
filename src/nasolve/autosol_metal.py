"""Metal-ion search by MR-SAD AutoSol directly after Phaser.

``nasolve autosol RUN --metal AG`` runs phenix.autosol in MR-SAD mode with the
original Phaser solution as the partial model, asking it to find sites of the
named metal(s). This is the step a crystallographer would otherwise run by
hand to see whether a metal ion is incorporated (for example in a
metal-mediated base pair) before the model is prepared further.

Contract
--------
* Requires a NASolve AutoMR run whose Phaser result is accepted
  (TFZ >= 8; 7 <= TFZ < 8 only with ``allow_mr_review``). It may run before
  or after PostMR: it always uses the original Phaser model, never a modified
  one, so a modeled ion cannot validate itself.
* The first metal is ``autosol.atom_type``; further metals go to
  ``autosol.phasing.mad_ha_add_list``. The partial-model RMSD
  (``autosol.phasing.partpdb_rms``) defaults to 1.0 A. The number of sites is
  left to AutoSol; model building is switched off.
* Each search writes a new ``AutoSolMetal/attempt_NNN`` directory; earlier
  attempts are never overwritten.
* The run's ``stage``/``status`` are unchanged. Results are added under
  ``report.json["metal_autosol"]`` with portable artifact references and
  SHA-256 checksums, and each attempt also has its own ``report.json``.
* Outcomes: ``METAL_SITES_FOUND`` (at least one site of a requested metal),
  ``METAL_NOT_FOUND`` (AutoSol finished without such a site - this is a normal
  result and usually means the ion is not incorporated), and
  ``METAL_AUTOSOL_FAILED`` (AutoSol stopped with an error or its outputs were
  unusable; often also just "no sites", but the log should be checked).
  None of these blocks the ordinary PostMR/AutoRefine path, and found sites
  are not added to the model automatically: they are for inspection in Coot.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

from .autosol import (
    AutoSolPreparationError,
    _effective_settings,
    _original_phaser_model,
    _pdb_atoms,
    _required_path,
    _run_mtz_dump,
    _sequence_source,
    _unit_cell,
    _verify_frozen_checksum,
    anomalous_labels,
    nearest_symmetry_distance,
    read_wavelength,
)
from .model_assessment import file_sha256
from .run_context import artifact_reference, resolve_artifact_path


DEFAULT_PARTIAL_MODEL_RMSD = 1.0
METAL_DIRECTORY = "AutoSolMetal"
SCHEMA_VERSION = 1

FOUND = "METAL_SITES_FOUND"
NOT_FOUND = "METAL_NOT_FOUND"
FAILED = "METAL_AUTOSOL_FAILED"

# Metals AutoSol can be asked to locate. Symbols are written in AutoSol's
# element casing. Non-metals (I, Br, Se, S, P) belong to the nucleotide
# heavy-atom branch, not to this search.
_METALS = (
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Th", "U",
)
METAL_SYMBOLS = {symbol.upper(): symbol for symbol in _METALS}

_NEAREST_RESIDUES = 3


@dataclass(frozen=True)
class MetalSearchResult:
    status: str
    message: str
    run_directory: Path
    attempt_directory: Path
    report_path: Path
    log_path: Path
    metals: tuple[str, ...]
    partial_model_rmsd: float
    phaser_model: Path
    phaser_map: Path
    heavy_atom_model: Path | None
    sites: tuple[dict[str, object], ...]


def parse_metals(values: Sequence[str]) -> tuple[str, ...]:
    """Normalize ``["ag", "Hg,Zn"]`` to ``("Ag", "Hg", "Zn")``.

    Order is kept (the first metal is AutoSol's primary atom type) and
    duplicates are dropped. Anything that is not a metal is rejected.
    """
    metals: list[str] = []
    for value in values:
        for item in str(value).replace(";", ",").replace(" ", ",").split(","):
            item = item.strip()
            if not item:
                continue
            symbol = METAL_SYMBOLS.get(item.upper())
            if symbol is None:
                raise AutoSolPreparationError(
                    f"{item!r} is not a metal AutoSol can search for "
                    "(nucleotide I/Br/Se use 'nasolve autosol' after PostMR)"
                )
            if symbol not in metals:
                metals.append(symbol)
    if not metals:
        raise AutoSolPreparationError("Name at least one metal, for example --metal AG")
    return tuple(metals)


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _read_report(run_directory: Path) -> tuple[Path, dict[str, object]]:
    selected = run_directory.expanduser().resolve()
    report_path = selected if selected.name == "report.json" else selected / "report.json"
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AutoSolPreparationError(f"Could not read NASolve report {report_path}: {exc}") from exc
    if not isinstance(report, dict):
        raise AutoSolPreparationError(f"NASolve report is not a JSON object: {report_path}")
    if report.get("workflow") != "automr":
        raise AutoSolPreparationError("The metal search requires a NASolve AutoMR run")
    return report_path.parent, report


def _phaser_record(report: Mapping[str, object]) -> Mapping[str, object]:
    execution = report.get("execution")
    phaser = execution.get("phaser") if isinstance(execution, Mapping) else None
    if not isinstance(phaser, Mapping):
        raise AutoSolPreparationError("This run has no Phaser result; run 'nasolve automr --execute' first")
    return phaser


def _check_phaser_accepted(report: Mapping[str, object], allow_mr_review: bool) -> float:
    """Return the TFZ, refusing failed MR and unapproved review results."""
    tfz = _phaser_record(report).get("best_tfz")
    if not isinstance(tfz, (int, float)) or isinstance(tfz, bool) or not math.isfinite(tfz):
        raise AutoSolPreparationError("The Phaser result has no TFZ; MR did not succeed")
    if tfz < 7.0:
        raise AutoSolPreparationError(f"Phaser TFZ {tfz:.2f} is below 7.0; MR did not succeed")
    if tfz < 8.0 and not allow_mr_review:
        raise AutoSolPreparationError(
            f"Phaser TFZ {tfz:.2f} needs review; pass --allow-mr-review to search anyway"
        )
    return float(tfz)


def _phaser_map(run: Path, report: Mapping[str, object]) -> Path:
    value = _phaser_record(report).get("solution_mtz")
    if value is not None:
        path = resolve_artifact_path(value, run)
        if path is not None and path.is_file():
            return path
        raise AutoSolPreparationError("The Phaser map (solution MTZ) recorded for this run is missing")
    fallback = run / "Phaser" / "mr_solution.mtz"
    if fallback.is_file():
        return fallback.resolve()
    raise AutoSolPreparationError("The Phaser map (solution MTZ) is missing")


def _next_attempt(run: Path) -> Path:
    root = run / METAL_DIRECTORY
    root.mkdir(exist_ok=True)
    numbers = []
    for child in root.iterdir():
        name = child.name
        if name.startswith("attempt_") and name[8:].isdigit():
            numbers.append(int(name[8:]))
    attempt = root / f"attempt_{(max(numbers) + 1 if numbers else 1):03d}"
    attempt.mkdir()  # never reuse an attempt directory
    return attempt


def _model_atoms(path: Path) -> list[dict[str, object]]:
    """Non-hydrogen atoms of the Phaser model with residue identity."""
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line[0:6] not in ("ATOM  ", "HETATM") or len(line) < 54:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        name = line[12:16].strip()
        if element in ("H", "D") or (not element and name[:1] in ("H", "D")):
            continue
        try:
            xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
        except ValueError:
            continue
        atoms.append({
            "chain": line[21:22].strip(),
            "residue_number": line[22:27].strip(),
            "residue_name": line[17:20].strip(),
            "atom_name": name,
            "coordinates": xyz,
        })
    return atoms


def _nearest_residues(
    site_xyz: Sequence[float],
    atoms: Sequence[Mapping[str, object]],
    cell: Sequence[float] | None,
    space_group: str | None,
) -> tuple[list[dict[str, object]], str]:
    """The closest atom of each of the nearest residues to one site."""
    symmetry = "crystallographic symmetry and lattice translations applied"
    use_symmetry = cell is not None and space_group is not None
    if use_symmetry:
        try:
            nearest_symmetry_distance(atoms[0]["coordinates"], site_xyz, cell, space_group)  # type: ignore[arg-type]
        except (AutoSolPreparationError, IndexError):
            use_symmetry = False
    if not use_symmetry:
        symmetry = (
            f"same asymmetric unit only (symmetry for {space_group!r} is not configured)"
            if space_group else "same asymmetric unit only (no CRYST1 record)"
        )
    best: dict[tuple[str, str, str], dict[str, object]] = {}
    for atom in atoms:
        xyz = atom["coordinates"]
        if use_symmetry:
            distance = float(nearest_symmetry_distance(xyz, site_xyz, cell, space_group)["distance"])  # type: ignore[arg-type]
        else:
            distance = math.dist(xyz, site_xyz)  # type: ignore[arg-type]
        key = (str(atom["chain"]), str(atom["residue_number"]), str(atom["residue_name"]))
        current = best.get(key)
        if current is None or distance < float(current["distance"]):
            best[key] = {
                "chain": key[0],
                "residue_number": key[1],
                "residue_name": key[2],
                "atom_name": atom["atom_name"],
                "distance": round(distance, 3),
            }
    ranked = sorted(best.values(), key=lambda item: float(item["distance"]))
    return ranked[:_NEAREST_RESIDUES], symmetry


def _heavy_atom_models(directory: Path) -> list[Path]:
    return sorted(path for path in directory.rglob("overall_best_ha_pdb.pdb") if path.is_file())


def _reference(path: Path, run: Path) -> dict[str, str]:
    reference = artifact_reference(path, run)
    reference["sha256"] = file_sha256(path)
    return reference


def _record(
    run: Path,
    report: dict[str, object],
    attempt: Path,
    payload: dict[str, object],
) -> Path:
    attempt_report = attempt / "report.json"
    _write_json(attempt_report, payload)
    block = report.get("metal_autosol")
    if not isinstance(block, dict):
        block = {"schema_version": SCHEMA_VERSION, "attempts": []}
    attempts = block.get("attempts")
    if not isinstance(attempts, list):
        attempts = []
    attempts.append(payload)
    block["attempts"] = attempts
    block["latest"] = payload["attempt"]
    report["metal_autosol"] = block
    _write_json(run / "report.json", report)
    return attempt_report


def execute_metal_search(
    run_directory: Path,
    autosol_executable: Path,
    mtz_dump_executable: Path,
    metals: Sequence[str],
    *,
    partial_model_rmsd: float = DEFAULT_PARTIAL_MODEL_RMSD,
    allow_mr_review: bool = False,
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> MetalSearchResult:
    """Search for metal sites with MR-SAD AutoSol guided by the Phaser model."""
    elements = parse_metals(metals)
    if not math.isfinite(partial_model_rmsd) or partial_model_rmsd <= 0:
        raise AutoSolPreparationError("The partial-model RMSD must be a positive number of Angstrom")

    # Everything that can be checked is checked before an attempt directory
    # exists, so a refused search leaves no trace in the run.
    run, report = _read_report(run_directory)
    tfz = _check_phaser_accepted(report, allow_mr_review)
    reflections = _required_path(report, "inputs", "reflections", run)
    summary = _required_path(report, "inputs", "summary", run)
    inputs = report.get("inputs")
    assert isinstance(inputs, Mapping)
    _verify_frozen_checksum(reflections, inputs.get("reflections_sha256"), "reflections")
    phaser_model = _original_phaser_model(run, report)
    phaser_map = _phaser_map(run, report)
    wavelength, wavelength_source = read_wavelength(summary)
    labels = anomalous_labels(_run_mtz_dump(reflections, mtz_dump_executable, environment))
    sequence_source, generated_sequence = _sequence_source(report, run)
    nproc = max(1, processor_count if processor_count is not None else (os.cpu_count() or 1))

    attempt = _next_attempt(run)
    log_path = attempt / "autosol.console.log"
    sequence_path = attempt / "sequence_input.txt"
    if sequence_source is not None:
        shutil.copyfile(sequence_source, sequence_path)
        sequence_origin = "run-frozen sequence file"
    else:
        sequence_path.write_text(str(generated_sequence), encoding="utf-8")
        sequence_origin = "post_mr_plan.sequences"

    # Full PHIL paths: short names such as atom_type also exist in the
    # wavelength/deriv scopes and are ambiguous to newer Phenix parsers.
    command = [
        str(autosol_executable.expanduser().resolve()),
        f"autosol.data={reflections}",
        f"autosol.seq_file={sequence_path}",
        "autosol.labels=" + " ".join(labels),
        f"autosol.atom_type={elements[0]}",
        f"autosol.lambda={wavelength:.8g}",
        f"autosol.phasing.input_partpdb_file={phaser_model}",
        f"autosol.phasing.partpdb_rms={partial_model_rmsd:g}",
    ]
    if len(elements) > 1:
        command.append("autosol.phasing.mad_ha_add_list=" + " ".join(elements[1:]))
    command.extend([
        "autosol.model_building.build=False",
        "autosol.model_building.phase_improve_and_build=False",
        f"autosol.general.nproc={nproc}",
    ])

    created = datetime.now(timezone.utc).isoformat()
    base: dict[str, object] = {
        "attempt": attempt.name,
        "created_utc": created,
        "metals": list(elements),
        "partial_model_rmsd": partial_model_rmsd,
        "phaser_tfz": tfz,
        "wavelength": wavelength,
        "wavelength_source": wavelength_source,
        "input_labels": list(labels),
        "sequence_source": sequence_origin,
        "inputs": {
            "reflections": _reference(reflections, run),
            "phaser_model": _reference(phaser_model, run),
            "phaser_map": _reference(phaser_map, run),
            "sequence": _reference(sequence_path, run),
        },
        "command": command,
        "directory": artifact_reference(attempt, run),
        "model_changed": False,
    }

    def finish(status: str, message: str, heavy_atom_model: Path | None,
               sites: list[dict[str, object]], extra: dict[str, object]) -> MetalSearchResult:
        payload = {**base, "status": status, "message": message, **extra}
        payload["outputs"] = {
            "console_log": _reference(log_path, run) if log_path.is_file() else None,
            "heavy_atom_model": _reference(heavy_atom_model, run) if heavy_atom_model else None,
        }
        payload["sites"] = sites
        report_path = _record(run, report, attempt, payload)
        return MetalSearchResult(
            status=status, message=message, run_directory=run, attempt_directory=attempt,
            report_path=report_path, log_path=log_path, metals=elements,
            partial_model_rmsd=partial_model_rmsd, phaser_model=phaser_model,
            phaser_map=phaser_map, heavy_atom_model=heavy_atom_model, sites=tuple(sites),
        )

    try:
        completed = subprocess.run(
            command, cwd=attempt,
            env=dict(environment) if environment is not None else None,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
        )
    except OSError as exc:
        log_path.write_text(f"NASolve could not launch AutoSol: {exc}\n", encoding="utf-8")
        return finish(FAILED, f"Could not launch AutoSol: {exc}", None, [], {"return_code": None})
    log_path.write_text(completed.stdout or "", encoding="utf-8")
    if completed.returncode:
        return finish(
            FAILED,
            f"AutoSol stopped with exit status {completed.returncode}. This often just means "
            "it could not place a metal site, but check the log.",
            None, [], {"return_code": completed.returncode},
        )

    effective: dict[str, object] = {}
    eff_files = sorted(attempt.rglob("autosol.eff"))
    if len(eff_files) == 1:
        try:
            effective = _effective_settings(eff_files[0], nproc)
            rmsd_values = re.findall(
                r"^\s*partpdb_rms\s*=\s*(\S+)\s*$",
                eff_files[0].read_text(encoding="utf-8", errors="replace"), re.I | re.M,
            )
            numeric = []
            for value in rmsd_values:
                try:
                    numeric.append(float(value))
                except ValueError:
                    continue
            if numeric and not any(math.isclose(v, partial_model_rmsd) for v in numeric):
                raise AutoSolPreparationError(
                    f"partpdb_rms is {numeric[0]:g}, not the requested {partial_model_rmsd:g}"
                )
            effective["partpdb_rms"] = numeric[0] if numeric else "not reported"
            effective["path"] = _reference(eff_files[0], run)
        except AutoSolPreparationError as exc:
            return finish(FAILED, f"AutoSol settings were not applied as requested: {exc}",
                          None, [], {"return_code": 0})
    extra: dict[str, object] = {"return_code": 0, "effective_parameters": effective}

    models = _heavy_atom_models(attempt)
    if len(models) > 1:
        return finish(FAILED, "AutoSol wrote more than one overall_best_ha_pdb.pdb; "
                      "inspect the attempt directory", None, [], extra)
    if not models:
        return finish(NOT_FOUND, "AutoSol placed no heavy-atom site. This is a normal result "
                      "and usually means the metal is not incorporated.", None, [], extra)
    heavy_atom_model = models[0]
    try:
        all_sites = _pdb_atoms(heavy_atom_model)
    except AutoSolPreparationError:
        all_sites = []
    wanted = {element.upper() for element in elements}
    found = [site for site in all_sites if str(site.get("element", "")).upper() in wanted]
    extra["other_sites"] = len(all_sites) - len(found)
    if not found:
        return finish(
            NOT_FOUND,
            "AutoSol placed no site of the requested metal(s). This is a normal result and "
            "usually means the metal is not incorporated.",
            heavy_atom_model, [], extra,
        )

    atoms = _model_atoms(phaser_model)
    try:
        *cell, space_group = _unit_cell(phaser_model)
    except AutoSolPreparationError:
        cell, space_group = None, None
    sites = []
    symmetry_note = ""
    for index, site in enumerate(found, start=1):
        nearest, symmetry_note = (
            _nearest_residues(site["coordinates"], atoms, cell, space_group)  # type: ignore[arg-type]
            if atoms else ([], "the Phaser model has no atoms")
        )
        sites.append({
            "index": index,
            "element": str(site["element"]).upper(),
            "atom_name": site["atom_name"],
            "coordinates": [round(value, 3) for value in site["coordinates"]],  # type: ignore[union-attr]
            "occupancy": site["occupancy"],
            "nearest_residues": nearest,
        })
    extra["distance_method"] = symmetry_note
    counts = ", ".join(
        f"{sum(1 for s in sites if s['element'] == e.upper())} {e}"
        for e in elements if any(s["element"] == e.upper() for s in sites)
    )
    return finish(
        FOUND,
        f"AutoSol found {counts} site(s); inspect them in Coot against the Phaser model and map",
        heavy_atom_model, sites, extra,
    )


__all__ = [
    "DEFAULT_PARTIAL_MODEL_RMSD",
    "FAILED",
    "FOUND",
    "METAL_DIRECTORY",
    "METAL_SYMBOLS",
    "MetalSearchResult",
    "NOT_FOUND",
    "execute_metal_search",
    "parse_metals",
]
