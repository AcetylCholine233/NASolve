"""Authoritative space-group checks for standard AutoMR frame recipes."""

from __future__ import annotations

import html
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping


class SymmetryError(RuntimeError):
    """Raised when standard-frame symmetry is missing, inconsistent, or unsafe."""


@dataclass(frozen=True)
class SpaceGroupEvidence:
    mtz_symbol: str
    mtz_number: int | None
    mtz_matrix_symbol: str | None
    cif_symbols: tuple[str, ...]
    summary_symbol: str | None  # None when the dataset has no summary.html
    normalized_class: str


@dataclass(frozen=True)
class StandardSymmetryAssessment:
    evidence: SpaceGroupEvidence
    mr_copies: int
    red_flag: str | None
    final_output_patch: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


_CIF_TAGS = (
    "_symmetry.space_group_name_H-M",
    "_space_group.name_H-M_alt",
)


def _clean_symbol(value: str) -> str:
    return value.strip().strip("'\"").strip()


def normalized_space_group(value: str) -> str:
    """Normalize notation while deliberately treating H3 and R3 as equivalent."""
    compact = re.sub(r"[\s:_]+", "", _clean_symbol(value)).upper()
    if compact in {"H3", "R3", "R3H"}:
        return "H3/R3"
    if compact == "P1":
        return "P1"
    return compact


def read_cif_space_groups(path: Path) -> tuple[str, ...]:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise SymmetryError(f"Could not read crystallographic metadata {path}: {exc}") from exc
    values: list[str] = []
    tag_pattern = "|".join(re.escape(tag) for tag in _CIF_TAGS)
    pattern = re.compile(
        rf"^(?:{tag_pattern})\s+(?:'([^']+)'|\"([^\"]+)\"|(\S+))",
        re.I,
    )
    for line in lines:
        match = pattern.match(line.strip())
        if match:
            value = next(group for group in match.groups() if group is not None)
            values.append(_clean_symbol(value))
    if not values:
        raise SymmetryError(f"No explicit space-group tag was found in {path.name}")
    return tuple(values)


def read_final_summary_space_group(path: Path) -> str:
    """Return only the last authoritative ``Spacegroup name`` entry."""
    try:
        raw = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise SymmetryError(f"Could not read {path}: {exc}") from exc
    matches: list[str] = []
    for raw_line in raw.splitlines():
        line = html.unescape(re.sub(r"<[^>]*>", " ", raw_line)).strip()
        match = re.fullmatch(r"Spacegroup\s+name\s+([A-Za-z][A-Za-z0-9 ]*)", line, re.I)
        if match:
            matches.append(_clean_symbol(match.group(1)))
    if not matches:
        raise SymmetryError(
            f"No final 'Spacegroup name' entry was found near the bottom of {path.name}"
        )
    return matches[-1]


def read_mtz_space_group(
    mtz_path: Path,
    mtz_dump_executable: Path,
    environment: Mapping[str, str] | None = None,
) -> tuple[str, int | None, str | None]:
    """Read the operational MTZ header using the discovered Phenix installation."""
    command = [str(mtz_dump_executable), str(mtz_path)]
    try:
        result = subprocess.run(
            command,
            env=dict(environment) if environment is not None else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise SymmetryError(f"Could not inspect MTZ space group: {exc}") from exc
    if result.returncode:
        raise SymmetryError(
            f"phenix.mtz.dump failed ({result.returncode}) while inspecting {mtz_path.name}"
        )
    symbol_match = re.search(
        r"^Space group symbol from file:\s*(.+?)\s*$", result.stdout, re.I | re.M
    )
    number_match = re.search(
        r"^Space group number from file:\s*(\d+)\s*$", result.stdout, re.I | re.M
    )
    matrix_match = re.search(
        r"^Space group from matrices:\s*(.+?)\s*$", result.stdout, re.I | re.M
    )
    if not symbol_match:
        raise SymmetryError(
            f"phenix.mtz.dump did not report a space-group symbol for {mtz_path.name}"
        )
    symbol = _clean_symbol(symbol_match.group(1))
    number = int(number_match.group(1)) if number_match else None
    matrix_symbol = None
    if matrix_match:
        matrix_symbol = re.sub(
            r"\s*\(No\.\s*\d+\)\s*$", "", _clean_symbol(matrix_match.group(1)), flags=re.I
        )
    return symbol, number, matrix_symbol


def assess_standard_symmetry(
    mtz_path: Path,
    cif_path: Path,
    summary_path: Path | None,
    mtz_dump_executable: Path,
    environment: Mapping[str, str] | None = None,
    allow_p1_standard: bool = False,
) -> StandardSymmetryAssessment:
    """Cross-check authoritative symmetry sources and apply the standard-frame gate.

    The MTZ and every explicit ``Data_1*.cif`` entry are always compared; the
    final ``summary.html`` entry joins the comparison when that file exists.
    """
    mtz_symbol, mtz_number, matrix_symbol = read_mtz_space_group(
        mtz_path, mtz_dump_executable, environment
    )
    cif_symbols = read_cif_space_groups(cif_path)
    summary_symbol = (
        read_final_summary_space_group(summary_path) if summary_path is not None else None
    )
    normalized_values = {
        normalized_space_group(mtz_symbol),
        *(normalized_space_group(value) for value in cif_symbols),
    }
    if summary_symbol is not None:
        normalized_values.add(normalized_space_group(summary_symbol))
    if matrix_symbol:
        normalized_values.add(normalized_space_group(matrix_symbol))
    if len(normalized_values) != 1:
        raise SymmetryError(
            "Authoritative space-group sources disagree: "
            f"MTZ={mtz_symbol!r}, CIF={list(cif_symbols)!r}, "
            + (f"final summary={summary_symbol!r}" if summary_symbol is not None
               else "no summary.html")
        )
    normalized = next(iter(normalized_values))
    evidence = SpaceGroupEvidence(
        mtz_symbol=mtz_symbol,
        mtz_number=mtz_number,
        mtz_matrix_symbol=matrix_symbol,
        cif_symbols=cif_symbols,
        summary_symbol=summary_symbol,
        normalized_class=normalized,
    )
    patch = {
        "required": True,
        "target_symbol": mtz_symbol,
        "reason": "restore the input space-group setting label after structure solution",
    }
    if normalized == "H3/R3":
        if mtz_number not in {None, 146}:
            raise SymmetryError(
                f"MTZ symbol {mtz_symbol!r} is H3/R3-equivalent but reports number {mtz_number}"
            )
        return StandardSymmetryAssessment(evidence, 1, None, patch)
    if normalized == "P1":
        if mtz_number not in {None, 1}:
            raise SymmetryError(
                f"MTZ symbol {mtz_symbol!r} is P1 but reports number {mtz_number}"
            )
        if not allow_p1_standard:
            raise SymmetryError(
                "Standard W/3GBI recipes are intended for H3/R3 data, but the "
                "authoritative inputs report P1. This path requires three MR copies "
                "and is strongly discouraged. Use --allow-p1-standard only after review."
            )
        return StandardSymmetryAssessment(
            evidence,
            3,
            "P1 standard-frame shunt enabled: Phaser must search for three copies",
            patch,
        )
    raise SymmetryError(
        f"Standard W/3GBI recipes require H3/R3 data; inputs report {mtz_symbol}"
    )
