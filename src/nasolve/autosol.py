"""Guarded MR-SAD AutoSol execution and anomalous-site validation."""

from __future__ import annotations

import html
import itertools
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Mapping, Sequence

from .model_assessment import file_sha256
from .reflection_data import ReflectionDataError, mtz_wavelength
from .run_context import artifact_reference, resolve_artifact_path


class AutoSolPreparationError(RuntimeError):
    """Raised when the conditional AutoSol branch cannot be run safely."""


@dataclass(frozen=True)
class AutoSolResult:
    status: str
    message: str
    run_directory: Path
    autosol_directory: Path
    report_path: Path
    log_path: Path
    heavy_atom_model: Path | None
    refinement_data: Path | None
    matched_distance: float | None


DEFAULT_SITE_MATCH_TOLERANCE = 4.0
DEFAULT_SITE_REVIEW_TOLERANCE = 8.0
_ATOM_TYPE = {"I": "I", "BR": "Br", "SE": "Se"}
_ANOMALOUS_LABEL_SETS = (
    ("I(+)", "SIGI(+)", "I(-)", "SIGI(-)"),
    ("F(+)", "SIGF(+)", "F(-)", "SIGF(-)"),
)


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _read_run_report(run_directory: Path) -> tuple[Path, dict[str, object]]:
    selected = run_directory.expanduser().resolve()
    report_path = selected if selected.name == "report.json" else selected / "report.json"
    run = report_path.parent
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AutoSolPreparationError(f"Could not read NASolve report {report_path}: {exc}") from exc
    if not isinstance(report, dict):
        raise AutoSolPreparationError(
            f"NASolve report is not a JSON object: {report_path}"
        )
    if report.get("workflow") != "automr":
        raise AutoSolPreparationError("AutoSol requires a NASolve AutoMR run")
    if report.get("stage") == "autosol" and report.get("status") == "AUTOSOL_READY":
        raise AutoSolPreparationError(
            f"AutoSol is already complete; refusing to overwrite {run / 'AutoSol'}"
        )
    if report.get("stage") != "postmr" or report.get("status") != "POSTMR_READY":
        raise AutoSolPreparationError("AutoSol requires a completed POSTMR_READY run")
    return run, report


def read_wavelength(summary_path: Path, required: bool = True) -> tuple[float, str] | None:
    """Read one consistent X-ray wavelength, preferring the high-precision assignment.

    With ``required=False`` a summary that states no wavelength returns None;
    conflicting or implausible values still stop.
    """
    try:
        raw = summary_path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise AutoSolPreparationError(f"Could not read wavelength source {summary_path}: {exc}") from exc
    plain = html.unescape(re.sub(r"<[^>]*>", " ", raw))
    assigned = [
        float(value)
        for value in re.findall(
            r"\bwavelength\b[^\n=]{0,80}=\s*([0-9]+(?:\.[0-9]+)?)",
            plain,
            re.I,
        )
    ]
    displayed = [
        float(value)
        for value in re.findall(
            r"\bWavelength\b\s+([0-9]+(?:\.[0-9]+)?)\s*(?:A|Å)\b",
            plain,
            re.I,
        )
    ]
    values = assigned or displayed
    if not values and not required:
        return None
    if not values:
        raise AutoSolPreparationError(
            f"No wavelength was found in {summary_path.name}"
        )
    if any(not 0.1 <= value <= 5.0 for value in values):
        raise AutoSolPreparationError(
            f"Implausible wavelength value in {summary_path.name}: {values}"
        )
    reference = assigned[0] if assigned else displayed[0]
    comparison = assigned + displayed
    if any(abs(value - reference) > 5e-5 for value in comparison):
        raise AutoSolPreparationError(
            f"Conflicting wavelength values in {summary_path.name}: {comparison}"
        )
    source = "high-precision summary assignment" if assigned else "summary display"
    return reference, source


# MTZ headers store the wavelength as a 4-byte float; summary text may round.
_WAVELENGTH_TOLERANCE = 1e-4


def resolve_wavelength(
    mtz_dump_output: str,
    reflections: Path,
    summary_path: Path | None,
) -> tuple[float, dict[str, object]]:
    """Return the MTZ header wavelength; an optional summary.html must agree."""
    try:
        wavelength, dataset = mtz_wavelength(mtz_dump_output)
    except ReflectionDataError as exc:
        raise AutoSolPreparationError(str(exc)) from exc
    source: dict[str, object] = {
        "method": "MTZ header (phenix.mtz.dump)",
        "path": str(reflections),
        "mtz_dataset": dataset,
        "summary_cross_check": None,
    }
    if summary_path is not None:
        found = read_wavelength(summary_path, required=False)
        if found is None:
            source["summary_cross_check"] = {
                "path": str(summary_path), "wavelength": None,
                "method": "summary states no wavelength; not cross-checked",
            }
            return wavelength, source
        summary_value, how = found
        if abs(summary_value - wavelength) > _WAVELENGTH_TOLERANCE:
            raise AutoSolPreparationError(
                f"The MTZ header wavelength ({wavelength:g} A) and "
                f"{summary_path.name} ({summary_value:g} A) disagree"
            )
        source["summary_cross_check"] = {
            "path": str(summary_path), "wavelength": summary_value, "method": how,
        }
    return wavelength, source


def _run_mtz_dump(
    mtz_path: Path,
    executable: Path,
    environment: Mapping[str, str] | None,
) -> str:
    try:
        completed = subprocess.run(
            [str(executable), str(mtz_path)],
            env=dict(environment) if environment is not None else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise AutoSolPreparationError(f"Could not inspect MTZ {mtz_path}: {exc}") from exc
    if completed.returncode:
        raise AutoSolPreparationError(
            f"phenix.mtz.dump failed while inspecting {mtz_path.name}"
        )
    return completed.stdout


def _contains_label(output: str, label: str) -> bool:
    return re.search(rf"(?<!\S){re.escape(label)}(?!\S)", output, re.I) is not None


def anomalous_labels(mtz_dump_output: str) -> tuple[str, str, str, str]:
    for labels in _ANOMALOUS_LABEL_SETS:
        if all(_contains_label(mtz_dump_output, label) for label in labels):
            return labels
    raise AutoSolPreparationError(
        "Input MTZ has no complete anomalous intensity or amplitude quartet"
    )


def _sequence_source(
    report: Mapping[str, object], run: Path
) -> tuple[Path | None, str | None]:
    inputs = report.get("inputs")
    if isinstance(inputs, Mapping) and inputs.get("frame_sequence") is not None:
        sequence_value = inputs["frame_sequence"]
        sequence = resolve_artifact_path(
            sequence_value,
            run,
            verify_checksum=False,
        )
        if sequence is None:
            raise AutoSolPreparationError("Frozen frame sequence is missing")
        if isinstance(sequence_value, Mapping):
            if sequence_value.get("anchor") not in {"run", "dataset", "repository"}:
                raise AutoSolPreparationError(
                    "Frozen frame sequence reference is not portable"
                )
            expected = sequence_value.get("sha256")
            if expected is None:
                raise AutoSolPreparationError(
                    "Frozen frame sequence reference has no checksum"
                )
        else:
            expected = inputs.get("frame_sequence_sha256")
        _verify_frozen_checksum(sequence, expected, "frame sequence")
        return sequence, "run-frozen frame sequence"
    model_value = inputs.get("model") if isinstance(inputs, Mapping) else None
    if isinstance(model_value, str):
        posix = PurePosixPath(model_value)
        windows = PureWindowsPath(model_value)
        if windows.is_absolute() and not posix.is_absolute():
            flavors = (windows,)
        elif posix.is_absolute() and not windows.is_absolute():
            flavors = (posix,)
        elif "\\" in model_value and "/" not in model_value:
            flavors = (windows,)
        elif "/" in model_value and "\\" not in model_value:
            flavors = (posix,)
        else:
            flavors = (posix, windows)
        legacy_parents = {
            str(path.parent / "seq_base.txt")
            for path in flavors
            if len(path.parts) > 1 and str(path.parent) not in {"", "."}
        }
        legacy_matches = {
            path
            for value in legacy_parents
            if (path := resolve_artifact_path(value, run)) is not None
        }
        if len(legacy_matches) == 1:
            return legacy_matches.pop(), "legacy frame-adjacent seq_base.txt"
        if len(legacy_matches) > 1:
            raise AutoSolPreparationError(
                "Legacy frame sequence resolves to multiple local files"
            )
    frame = report.get("frame")
    catalogue = frame.get("catalogue_directory") if isinstance(frame, Mapping) else None
    if isinstance(catalogue, str) and catalogue:
        if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", catalogue) is None:
            raise AutoSolPreparationError(
                "Run report contains an unsafe frame catalogue directory"
            )
        sequence = resolve_artifact_path(
            {
                "anchor": "repository",
                "relative_path": f"MR_frames/{catalogue}/seq_base.txt",
            },
            run,
        )
        if sequence is not None:
            return sequence, "declared frame seq_base.txt"
    if model_value is not None:
        model = resolve_artifact_path(model_value, run)
        adjacent = model.parent / "seq_base.txt" if model is not None else None
        if adjacent is not None and adjacent.is_file():
            return adjacent, "frame-adjacent seq_base.txt"
    plan = report.get("post_mr_plan")
    if isinstance(plan, Mapping):
        sequences = plan.get("sequences")
        if isinstance(sequences, Mapping) and sequences and all(
            isinstance(value, str) and value for value in sequences.values()
        ):
            text = "\n\n".join(str(value) for value in sequences.values()) + "\n"
            return None, text
    raise AutoSolPreparationError(
        "AutoSol needs MR_frames/FRAME/seq_base.txt or complete run sequences"
    )


def _original_phaser_model(run: Path, report: Mapping[str, object]) -> Path:
    execution = report.get("execution")
    if isinstance(execution, Mapping):
        phaser = execution.get("phaser")
        if isinstance(phaser, Mapping) and isinstance(phaser.get("solution_pdb"), str):
            model = resolve_artifact_path(phaser["solution_pdb"], run)
            if model is not None:
                return model
    fallback = run / "Phaser" / "mr_solution.pdb"
    if fallback.is_file():
        return fallback.resolve()
    raise AutoSolPreparationError("The original Phaser MR model is missing")


def _required_path(
    report: Mapping[str, object], section: str, key: str, run: Path
) -> Path:
    selected = report.get(section)
    value = selected.get(key) if isinstance(selected, Mapping) else None
    if not isinstance(value, str):
        raise AutoSolPreparationError(f"Run report has no {section}.{key} path")
    path = resolve_artifact_path(value, run)
    if path is None:
        raise AutoSolPreparationError(f"Required run input is missing: {value}")
    return path


def _optional_path(
    report: Mapping[str, object], section: str, key: str, run: Path
) -> Path | None:
    """A recorded optional input (None when the dataset did not have it)."""
    selected = report.get(section)
    if isinstance(selected, Mapping) and selected.get(key) is None:
        return None
    return _required_path(report, section, key, run)


def _verify_frozen_checksum(path: Path, expected: object, description: str) -> None:
    if expected is None:
        return
    if not isinstance(expected, str) or re.fullmatch(r"[0-9a-fA-F]{64}", expected) is None:
        raise AutoSolPreparationError(f"Frozen {description} checksum is malformed")
    if file_sha256(path) != expected.casefold():
        raise AutoSolPreparationError(
            f"Frozen {description} does not match its recorded checksum"
        )


def _effective_settings(path: Path, expected_nproc: int) -> dict[str, object]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise AutoSolPreparationError(f"Could not inspect AutoSol settings {path}: {exc}") from exc
    if not re.search(r"^\s*build\s*=\s*False\s*$", text, re.I | re.M):
        raise AutoSolPreparationError("AutoSol effective settings do not guarantee build=False")
    if not re.search(
        r"^\s*phase_improve_and_build\s*=\s*False\s*$", text, re.I | re.M
    ):
        raise AutoSolPreparationError(
            "AutoSol effective settings do not guarantee phase_improve_and_build=False"
        )
    site_values = re.findall(r"^\s*sites\s*=\s*(\S+)\s*$", text, re.I | re.M)
    if not site_values or any(value.casefold() not in {"none", "null"} for value in site_values):
        raise AutoSolPreparationError("AutoSol effective settings contain specified HA sites")
    nproc_values = [
        int(value) for value in re.findall(r"^\s*nproc\s*=\s*(\d+)\s*$", text, re.I | re.M)
    ]
    if expected_nproc not in nproc_values:
        raise AutoSolPreparationError(
            f"AutoSol effective settings do not contain nproc={expected_nproc}"
        )
    return {
        "build": False,
        "phase_improve_and_build": False,
        "sites": None,
        "nproc": expected_nproc,
    }


def _find_unique(root: Path, name: str) -> Path:
    candidates = [
        path for path in root.rglob(name)
        if path.is_file() and "PDS" not in path.parts and not path.name.endswith(".str")
    ]
    if len(candidates) != 1:
        raise AutoSolPreparationError(
            f"Expected one {name} beneath {root}, found {len(candidates)}"
        )
    return candidates[0].resolve()


def _hl_labels(output: str) -> bool:
    normalized = output.upper().replace("HL_A", "HLA").replace("HL_B", "HLB")
    normalized = normalized.replace("HL_C", "HLC").replace("HL_D", "HLD")
    # Phenix 1.20.1 writes model-combined coefficients as
    # HLAM/HLBM/HLCM/HLDM; newer files often use HLA/HLB/HLC/HLD.
    return all(
        re.search(rf"\b{label}(?:M)?\b", normalized)
        for label in ("HLA", "HLB", "HLC", "HLD")
    )


def _discover_refinement_data(
    root: Path,
    mtz_dump_executable: Path,
    environment: Mapping[str, str] | None,
) -> tuple[Path, list[str]]:
    candidates = [
        path for path in root.rglob("*.mtz")
        if path.is_file() and "PDS" not in path.parts and "TEMP" not in path.parts
    ]
    candidates.sort(key=lambda path: (
        path.name != "overall_best_refine_data.mtz",
        "with_hl_anom" not in path.name,
        len(path.parts),
        str(path),
    ))
    inspected: list[str] = []
    for candidate in candidates:
        output = _run_mtz_dump(candidate, mtz_dump_executable, environment)
        if not _hl_labels(output):
            continue
        try:
            labels = anomalous_labels(output)
        except AutoSolPreparationError:
            continue
        inspected.append(str(candidate.resolve()))
        return candidate.resolve(), list(labels)
    raise AutoSolPreparationError(
        "AutoSol produced no MTZ containing both Hendrickson-Lattman phases "
        "and anomalous observations"
    )


def _density_modified_map_outputs(root: Path, run: Path) -> dict[str, object]:
    """Record optional view-map provenance without changing phase acceptance."""
    outputs: dict[str, object] = {
        "density_modified_map": None,
        "density_modified_map_sha256": None,
        "density_modified_map_status": "missing",
        "density_modified_map_diagnostic": "No density-modified map was produced",
    }
    try:
        for pattern in (
            "overall_best_denmod_map_coeffs.mtz", "resolve_histograms*.mtz"
        ):
            candidates = sorted(
                path for path in root.rglob(pattern)
                if path.is_file()
                and not {"PDS", "TEMP"}.intersection(path.relative_to(root).parts)
            )
            if not candidates:
                continue
            if len(candidates) != 1:
                outputs["density_modified_map_status"] = "ambiguous"
                outputs["density_modified_map_diagnostic"] = (
                    f"Expected one density-modified map matching {pattern}; "
                    f"found {len(candidates)}: "
                    + ", ".join(path.relative_to(root).as_posix() for path in candidates)
                )
                return outputs
            selected = candidates[0].resolve()
            if not selected.is_relative_to(root.resolve()):
                raise ValueError("density-modified map is outside the AutoSol output directory")
            checksum = file_sha256(selected)
            outputs.update({
                "density_modified_map": {
                    **artifact_reference(selected, run), "sha256": checksum,
                },
                "density_modified_map_sha256": checksum,
                "density_modified_map_status": "available",
                "density_modified_map_diagnostic": None,
            })
            return outputs
    except (OSError, ValueError) as exc:
        outputs["density_modified_map_status"] = "unavailable"
        outputs["density_modified_map_diagnostic"] = (
            f"Could not record the density-modified map: {exc}"
        )
    return outputs


def _pdb_atoms(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise AutoSolPreparationError(f"Could not read HA model {path}: {exc}") from exc
    for line in lines:
        if line[0:6].strip().upper() not in {"ATOM", "HETATM"} or len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            letters = "".join(character for character in atom_name if character.isalpha()).upper()
            element = next((item for item in ("BR", "SE") if letters.startswith(item)), letters[:1])
        try:
            coordinates = tuple(float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
            occupancy = float(line[54:60]) if len(line) >= 60 and line[54:60].strip() else 1.0
        except ValueError as exc:
            raise AutoSolPreparationError(f"Invalid HA record in {path.name}: {line}") from exc
        atoms.append({
            "atom_name": atom_name,
            "element": element,
            "coordinates": coordinates,
            "occupancy": occupancy,
        })
    if not atoms:
        raise AutoSolPreparationError(f"AutoSol HA model contains no atoms: {path}")
    return atoms


def _unit_cell(model_path: Path) -> tuple[float, float, float, float, float, float, str]:
    try:
        lines = model_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise AutoSolPreparationError(f"Could not read model unit cell: {exc}") from exc
    line = next((item for item in lines if item.startswith("CRYST1")), None)
    if line is None:
        raise AutoSolPreparationError("Prepared model has no CRYST1 record for HA-site matching")
    try:
        values = tuple(float(line[start:end]) for start, end in (
            (6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54)
        ))
    except ValueError as exc:
        raise AutoSolPreparationError("Prepared model has an invalid CRYST1 record") from exc
    return (*values, line[55:66].strip())


def _orthogonalization_matrix(
    cell: Sequence[float],
) -> tuple[tuple[float, float, float], ...]:
    a, b, c, alpha_deg, beta_deg, gamma_deg = cell
    alpha, beta, gamma = map(math.radians, (alpha_deg, beta_deg, gamma_deg))
    sin_gamma = math.sin(gamma)
    if abs(sin_gamma) < 1e-8:
        raise AutoSolPreparationError("Invalid unit-cell gamma angle")
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sin_gamma
    cz_squared = c * c - cx * cx - cy * cy
    if cz_squared <= 0:
        raise AutoSolPreparationError("Invalid unit-cell geometry")
    return (
        (a, b * math.cos(gamma), cx),
        (0.0, b * sin_gamma, cy),
        (0.0, 0.0, math.sqrt(cz_squared)),
    )


def _matrix_vector(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> tuple[float, float, float]:
    return tuple(sum(row[index] * vector[index] for index in range(3)) for row in matrix)  # type: ignore[return-value]


def _inverse_upper(matrix: Sequence[Sequence[float]]) -> tuple[tuple[float, float, float], ...]:
    m00, m01, m02 = matrix[0]
    _, m11, m12 = matrix[1]
    _, _, m22 = matrix[2]
    return (
        (1.0 / m00, -m01 / (m00 * m11), (m01 * m12 - m02 * m11) / (m00 * m11 * m22)),
        (0.0, 1.0 / m11, -m12 / (m11 * m22)),
        (0.0, 0.0, 1.0 / m22),
    )


def _symmetry_operators(space_group: str) -> tuple[tuple[str, object], ...]:
    compact = re.sub(r"[\s:_]+", "", space_group).upper()
    if compact in {"H3", "R3", "R3H"}:
        return (
            ("x,y,z", lambda x, y, z: (x, y, z)),
            ("-y,x-y,z", lambda x, y, z: (-y, x - y, z)),
            ("-x+y,-x,z", lambda x, y, z: (-x + y, -x, z)),
        )
    if compact == "P1":
        return (("x,y,z", lambda x, y, z: (x, y, z)),)
    raise AutoSolPreparationError(
        f"Symmetry-aware HA matching is not configured for space group {space_group!r}"
    )


def nearest_symmetry_distance(
    model_xyz: Sequence[float],
    site_xyz: Sequence[float],
    cell: Sequence[float],
    space_group: str,
) -> dict[str, object]:
    """Return the closest lattice/symmetry equivalent of ``model_xyz``."""
    orth = _orthogonalization_matrix(cell)
    fractional = _matrix_vector(_inverse_upper(orth), model_xyz)
    site_fractional = _matrix_vector(_inverse_upper(orth), site_xyz)
    compact = re.sub(r"[\s:_]+", "", space_group).upper()
    centering = (
        ((0.0, 0.0, 0.0), (2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0))
        if compact in {"H3", "R3", "R3H"}
        else ((0.0, 0.0, 0.0),)
    )
    best: dict[str, object] | None = None
    for operator_name, operator in _symmetry_operators(space_group):
        transformed = operator(*fractional)  # type: ignore[operator]
        for center in centering:
            centered = tuple(transformed[index] + center[index] for index in range(3))
            approximate = tuple(round(site_fractional[index] - centered[index]) for index in range(3))
            for offset in itertools.product((-1, 0, 1), repeat=3):
                lattice = tuple(int(approximate[index] + offset[index]) for index in range(3))
                equivalent_fractional = tuple(centered[index] + lattice[index] for index in range(3))
                equivalent_xyz = _matrix_vector(orth, equivalent_fractional)
                distance = math.dist(equivalent_xyz, site_xyz)
                if best is None or distance < float(best["distance"]):
                    best = {
                        "distance": distance,
                        "operator": operator_name,
                        "centering_translation": list(center),
                        "lattice_translation": list(lattice),
                        "equivalent_model_coordinates": list(equivalent_xyz),
                    }
    if best is None:
        raise AutoSolPreparationError("No crystallographic operators were available")
    return best


def _match_sites(
    candidates: Sequence[Mapping[str, object]],
    sites: Sequence[Mapping[str, object]],
    model_path: Path,
    automatic_tolerance: float,
    review_tolerance: float,
) -> list[dict[str, object]]:
    a, b, c, alpha, beta, gamma, space_group = _unit_cell(model_path)
    cell = (a, b, c, alpha, beta, gamma)
    attempts: list[dict[str, object]] = []
    for candidate in candidates:
        model_xyz = candidate.get("coordinates")
        element = str(candidate.get("element", "")).upper()
        if not isinstance(model_xyz, list) or len(model_xyz) != 3:
            raise AutoSolPreparationError("PostMR anomalous candidate has invalid coordinates")
        compatible = [site for site in sites if str(site.get("element", "")).upper() == element]
        if not compatible:
            attempts.append({"candidate": dict(candidate), "matched": False, "reason": "no same-element HA site"})
            continue
        matches: list[dict[str, object]] = []
        for site in compatible:
            site_xyz = site["coordinates"]
            match = nearest_symmetry_distance(model_xyz, site_xyz, cell, space_group)  # type: ignore[arg-type]
            matches.append({
                **match,
                "site_coordinates": list(site_xyz),  # type: ignore[arg-type]
                "site_occupancy": site["occupancy"],
            })
        nearest = min(matches, key=lambda item: float(item["distance"]))
        distance = float(nearest["distance"])
        result = {
            "candidate": dict(candidate),
            "within_automatic_limit": distance <= automatic_tolerance,
            "within_review_limit": distance <= review_tolerance,
            "nearest_site": {
                **nearest,
                "distance": round(distance, 3),
            },
        }
        attempts.append(result)
    return attempts


def _autosol_outcome(
    best_distance: float | None,
    automatic_tolerance: float,
    review_tolerance: float,
) -> tuple[str, str, bool]:
    if best_distance is not None and best_distance <= automatic_tolerance:
        return (
            "AUTOSOL_READY",
            "MR-SAD phases and a model-associated anomalous site are ready for refinement",
            True,
        )
    if best_distance is not None and best_distance <= review_tolerance:
        return (
            "AUTOSOL_REVIEW",
            f"Nearest anomalous site is {best_distance:.2f} A from the expected atom; "
            "continue, but inspect the site before using AutoSol phases",
            False,
        )
    detail = (
        f"nearest same-element site is {best_distance:.2f} A away"
        if best_distance is not None
        else "no same-element site could be matched"
    )
    return (
        "AUTOSOL_WARNING",
        f"AutoSol did not validate the expected nucleotide heavy atom ({detail}); "
        "continue without AutoSol phases and review the anomalous signal",
        False,
    )


def _record_result(
    run: Path,
    report: dict[str, object],
    autosol_directory: Path,
    payload: dict[str, object],
) -> Path:
    report_path = autosol_directory / "report.json"
    _write_json(report_path, payload)
    report["stage"] = "autosol"
    report["status"] = payload["status"]
    report["message"] = payload["message"]
    report["updated_utc"] = payload["created_utc"]
    report["autosol"] = payload
    _write_json(run / "report.json", report)
    return report_path


def _warning_result(
    run: Path,
    report: dict[str, object],
    autosol_directory: Path,
    reason: str,
    *,
    log_path: Path | None = None,
    command: list[str] | None = None,
) -> AutoSolResult:
    message = (
        "AutoSol did not produce validated experimental phases; continue with the "
        "MR refinement path and inspect this warning"
    )
    payload: dict[str, object] = {
        "status": "AUTOSOL_WARNING",
        "message": message,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "failure_reason": reason,
        "command": command,
        "outputs": {
            "heavy_atom_model": None,
            "refinement_data": None,
            "console_log": str(log_path) if log_path else None,
        },
        "use_for_refinement": False,
        "continuation": "continue without AutoSol phases; user review requested",
    }
    report_path = _record_result(run, report, autosol_directory, payload)
    return AutoSolResult(
        status="AUTOSOL_WARNING",
        message=message,
        run_directory=run,
        autosol_directory=autosol_directory,
        report_path=report_path,
        log_path=log_path or autosol_directory / "autosol.console.log",
        heavy_atom_model=None,
        refinement_data=None,
        matched_distance=None,
    )


def execute_autosol(
    run_directory: Path,
    autosol_executable: Path,
    mtz_dump_executable: Path,
    *,
    environment: Mapping[str, str] | None = None,
    match_tolerance: float = DEFAULT_SITE_MATCH_TOLERANCE,
    review_tolerance: float = DEFAULT_SITE_REVIEW_TOLERANCE,
    processor_count: int | None = None,
) -> AutoSolResult:
    """Run the conditional MR-SAD branch and validate a model-associated HA site."""
    run, report = _read_run_report(run_directory)
    postmr = report.get("postmr")
    anomalous = postmr.get("anomalous") if isinstance(postmr, Mapping) else None
    candidates = anomalous.get("candidates") if isinstance(anomalous, Mapping) else None
    if not isinstance(candidates, list) or not candidates:
        raise AutoSolPreparationError(
            "PostMR found no nucleotide heavy atom; AutoSol is not required"
        )
    if not all(isinstance(candidate, Mapping) for candidate in candidates):
        raise AutoSolPreparationError("PostMR anomalous candidates are malformed")
    elements = {str(candidate.get("element", "")).upper() for candidate in candidates}
    if len(elements) != 1 or next(iter(elements)) not in _ATOM_TYPE:
        raise AutoSolPreparationError(
            "One AutoSol run currently requires exactly one supported HA element; "
            f"found {sorted(elements)}"
        )
    element = next(iter(elements))
    autosol_directory = run / "AutoSol"
    try:
        autosol_directory.mkdir()
    except FileExistsError as exc:
        raise AutoSolPreparationError(
            f"AutoSol directory already exists; refusing to overwrite {autosol_directory}"
        ) from exc
    log_path = autosol_directory / "autosol.console.log"
    try:
        reflections = _required_path(report, "inputs", "reflections", run)
        summary = _optional_path(report, "inputs", "summary", run)
        inputs = report.get("inputs")
        assert isinstance(inputs, Mapping)
        _verify_frozen_checksum(
            reflections, inputs.get("reflections_sha256"), "reflections"
        )
        prepared_value = postmr.get("prepared_model") if isinstance(postmr, Mapping) else None
        prepared_model = resolve_artifact_path(prepared_value, run)
        if prepared_model is None:
            raise AutoSolPreparationError("PostMR prepared model is missing")
        _verify_frozen_checksum(
            prepared_model,
            postmr.get("prepared_sha256") if isinstance(postmr, Mapping) else None,
            "PostMR prepared model",
        )
        phaser_model = _original_phaser_model(run, report)
        input_dump = _run_mtz_dump(reflections, mtz_dump_executable, environment)
        wavelength, wavelength_source = resolve_wavelength(input_dump, reflections, summary)
        labels = anomalous_labels(input_dump)
        sequence_source, generated_sequence = _sequence_source(report, run)
        nproc = max(1, processor_count if processor_count is not None else (os.cpu_count() or 1))
    except AutoSolPreparationError as exc:
        return _warning_result(
            run, report, autosol_directory, str(exc), log_path=log_path
        )
    sequence_path = autosol_directory / "sequence_input.txt"
    if sequence_source is not None:
        shutil.copyfile(sequence_source, sequence_path)
        sequence_origin = str(sequence_source)
    else:
        sequence_path.write_text(str(generated_sequence), encoding="utf-8")
        sequence_origin = "post_mr_plan.sequences"

    # Use the full AutoSol PHIL paths. Leaf names such as data also exist in
    # wavelength/native/deriv scopes and are ambiguous in newer Phenix parsers.
    command = [
        str(autosol_executable.expanduser().resolve()),
        f"autosol.data={reflections}",
        f"autosol.seq_file={sequence_path}",
        "autosol.labels=" + " ".join(labels),
        f"autosol.atom_type={_ATOM_TYPE[element]}",
        f"autosol.lambda={wavelength:.8g}",
        f"autosol.phasing.input_partpdb_file={phaser_model}",
        "autosol.model_building.build=False",
        "autosol.model_building.phase_improve_and_build=False",
        f"autosol.general.nproc={nproc}",
    ]
    try:
        completed = subprocess.run(
            command,
            cwd=autosol_directory,
            env=dict(environment) if environment is not None else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except OSError as exc:
        log_path.write_text(f"NASolve could not launch AutoSol: {exc}\n", encoding="utf-8")
        return _warning_result(
            run,
            report,
            autosol_directory,
            f"Could not launch AutoSol: {exc}",
            log_path=log_path,
            command=command,
        )
    log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode:
        return _warning_result(
            run,
            report,
            autosol_directory,
            f"AutoSol exited with status {completed.returncode}",
            log_path=log_path,
            command=command,
        )
    try:
        effective_path = _find_unique(autosol_directory, "autosol.eff")
        effective = _effective_settings(effective_path, nproc)
        heavy_atom_model = _find_unique(autosol_directory, "overall_best_ha_pdb.pdb")
        refinement_data, phase_labels = _discover_refinement_data(
            autosol_directory, mtz_dump_executable, environment
        )
        sites = _pdb_atoms(heavy_atom_model)
        attempts = _match_sites(
            candidates,
            sites,
            prepared_model,
            match_tolerance,
            review_tolerance,
        )  # type: ignore[arg-type]
    except AutoSolPreparationError as exc:
        return _warning_result(
            run,
            report,
            autosol_directory,
            str(exc),
            log_path=log_path,
            command=command,
        )
    comparable = [
        attempt for attempt in attempts
        if isinstance(attempt.get("nearest_site"), Mapping)
    ]
    best = min(
        comparable,
        key=lambda item: float(item["nearest_site"]["distance"]),  # type: ignore[index]
    ) if comparable else None
    best_distance = (
        float(best["nearest_site"]["distance"])  # type: ignore[index]
        if best is not None else None
    )
    status, message, use_for_refinement = _autosol_outcome(
        best_distance, match_tolerance, review_tolerance
    )
    automatically_accepted = [
        attempt for attempt in attempts if attempt.get("within_automatic_limit") is True
    ]
    review_candidates = [
        attempt for attempt in attempts
        if attempt.get("within_review_limit") is True
        and attempt.get("within_automatic_limit") is not True
    ]
    payload = {
        "status": status,
        "message": message,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "element": element,
        "wavelength": wavelength,
        "wavelength_source": wavelength_source,
        "input_reflections": str(reflections),
        "input_labels": list(labels),
        "input_model": str(phaser_model),
        "input_model_role": "original Phaser MR model; not the iodine-containing PostMR model",
        "sequence": {
            "source": sequence_origin,
            "copied_input": str(sequence_path),
        },
        "command": command,
        "effective_parameters": {
            "path": str(effective_path),
            **effective,
            "site_count_argument_supplied": False,
        },
        "outputs": {
            "heavy_atom_model": str(heavy_atom_model),
            "heavy_atom_model_sha256": file_sha256(heavy_atom_model),
            "refinement_data": str(refinement_data),
            "refinement_data_sha256": file_sha256(refinement_data),
            "refinement_data_phase_labels": phase_labels,
            **_density_modified_map_outputs(autosol_directory, run),
            "console_log": str(log_path),
        },
        "site_validation": {
            "automatic_acceptance_distance": match_tolerance,
            "review_distance": review_tolerance,
            "acceptance_rule": (
                "<= automatic distance: use phases; <= review distance: request inspection; "
                "otherwise continue without phases"
            ),
            "attempts": attempts,
            "accepted": automatically_accepted,
            "review_candidates": review_candidates,
            "best_distance": best_distance,
        },
        "use_for_refinement": use_for_refinement,
        "continuation": (
            "use validated AutoSol refinement data"
            if use_for_refinement
            else "continue without AutoSol phases pending user review"
        ),
    }
    report_path = _record_result(run, report, autosol_directory, payload)
    return AutoSolResult(
        status=status,
        message=message,
        run_directory=run,
        autosol_directory=autosol_directory,
        report_path=report_path,
        log_path=log_path,
        heavy_atom_model=heavy_atom_model,
        refinement_data=refinement_data,
        matched_distance=best_distance,
    )


__all__ = [
    "AutoSolPreparationError",
    "AutoSolResult",
    "DEFAULT_SITE_MATCH_TOLERANCE",
    "DEFAULT_SITE_REVIEW_TOLERANCE",
    "anomalous_labels",
    "execute_autosol",
    "nearest_symmetry_distance",
    "read_wavelength",
    "resolve_wavelength",
]
