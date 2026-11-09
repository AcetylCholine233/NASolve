"""Choose which observed reflection array Phaser reads from the input MTZ.

``phenix.phaser`` started from the command line without ``phaser.labin``
chooses amplitudes whenever the file has any (``F,SIGF`` first) and only falls
back to intensities when it has none. Phaser's intensity likelihood target is
the better choice for weak or anisotropically truncated (STARANISO) data, and
it is what the Phenix GUI uses, so AutoMR can request intensities explicitly.

``phaser.labin`` must equal the label string that ``iotbx`` gives the array,
for example ``IMEAN,SIGIMEAN``. The labels are therefore chosen from the MTZ
column *types* reported by ``phenix.mtz.dump`` (``J`` intensity followed by
``Q`` standard deviation), following the way ``iotbx.mtz`` groups columns,
rather than from a guessed column name.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping


AUTO = "auto"
INTENSITIES = "intensities"
PHASER_DATA_CHOICES = (AUTO, INTENSITIES)
_ALIASES = {
    "auto": AUTO,
    "default": AUTO,
    "intensities": INTENSITIES,
    "intensity": INTENSITIES,
    "i": INTENSITIES,
    "isigi": INTENSITIES,
    "i/sigi": INTENSITIES,
    "i,sigi": INTENSITIES,
}
SCHEMA_VERSION = 1


class ReflectionDataError(RuntimeError):
    """Raised when the requested Phaser data cannot be selected safely."""


@dataclass(frozen=True)
class MtzColumn:
    dataset: int  # ordinal of the MTZ dataset in phenix.mtz.dump order
    label: str
    type: str


def normalize_phaser_data(value: object) -> str:
    """Return ``auto`` or ``intensities`` for a user setting."""
    if value is None:
        return AUTO
    if not isinstance(value, str):
        raise ReflectionDataError("phaser_data must be auto or intensities")
    key = "".join(value.split()).casefold()
    if key not in _ALIASES:
        raise ReflectionDataError(
            f"phaser_data must be auto or intensities (I), not {value!r}"
        )
    return _ALIASES[key]


# A column row of ``phenix.mtz.dump``'s summary, for example
# ``    IMEAN   1635 100.00%  -93.85  4926.23 J: intensity``
_COLUMN_RE = re.compile(
    r"^\s+(\S+)\s+\d+\s+[0-9.]+%\s+\S+\s+\S+\s+([A-Za-z]):(?:\s|$)"
)
_DATASET_RE = re.compile(r"^\s*Dataset\s+(\d+)\s*:\s*$")


def parse_mtz_dump_columns(text: str) -> tuple[MtzColumn, ...]:
    """Return the column labels and types listed by ``phenix.mtz.dump``."""
    columns: list[MtzColumn] = []
    # Dataset numbers restart in every crystal, so number them in file order.
    dataset = 0
    for line in text.splitlines():
        if _DATASET_RE.match(line):
            dataset += 1
            continue
        match = _COLUMN_RE.match(line)
        if match:
            columns.append(MtzColumn(dataset, match.group(1), match.group(2).upper()))
    return tuple(columns)


def _is_anomalous_label(label: str) -> bool:
    upper = label.upper()
    return (
        "(+)" in upper or "(-)" in upper
        or upper.endswith(("+", "-", "PLUS", "MINUS"))
    )


def mean_intensity_label_strings(columns: Iterable[MtzColumn]) -> list[str]:
    """Return ``iotbx`` label strings of merged-intensity (J,Q) arrays, in order.

    Columns are grouped per MTZ dataset in file order. ``F Q D Q`` starts a
    reconstructed-amplitude group and consumes all four columns, exactly as
    ``iotbx.mtz.object.group_columns`` does, so a J column is only paired with
    the Q column that ``iotbx`` itself pairs it with.
    """
    by_dataset: dict[int, list[MtzColumn]] = {}
    for column in columns:
        by_dataset.setdefault(column.dataset, []).append(column)
    found: list[str] = []
    for dataset_columns in by_dataset.values():
        types = "".join(column.type for column in dataset_columns)
        index = 0
        while index < len(dataset_columns):
            kind = types[index]
            if types[index:index + 4] == "FQDQ":
                index += 4
                continue
            if kind in "JFED" and index + 1 < len(types) and types[index + 1] in "QP":
                first, second = dataset_columns[index], dataset_columns[index + 1]
                if (
                    kind == "J" and types[index + 1] == "Q"
                    and not _is_anomalous_label(first.label)
                    and not _is_anomalous_label(second.label)
                ):
                    found.append(f"{first.label},{second.label}")
                index += 2
                continue
            if kind in "GK":
                if types[index:index + 4] in {"GLGL", "GPGP", "KMKM", "GGLL", "GGPP", "KKMM"}:
                    index += 4
                    continue
                if types[index:index + 2] in {"GG", "KK"}:
                    index += 2
                    continue
            index += 1
    return found


def select_intensity_labels(columns: Iterable[MtzColumn]) -> str:
    """Return the first merged intensity array, or explain why there is none."""
    columns = tuple(columns)
    if not columns:
        raise ReflectionDataError(
            "phenix.mtz.dump did not list any MTZ columns, so intensity columns "
            "could not be identified"
        )
    candidates = mean_intensity_label_strings(columns)
    if not candidates:
        present = ", ".join(f"{column.label}({column.type})" for column in columns)
        raise ReflectionDataError(
            "phaser_data = intensities needs a merged intensity column (MTZ type J) "
            "followed by its standard deviation (type Q); the MTZ has: " + present
        )
    return candidates[0]


def read_mtz_dump(
    mtz_path: Path,
    mtz_dump_executable: Path,
    environment: Mapping[str, str] | None = None,
) -> str:
    try:
        completed = subprocess.run(
            [str(mtz_dump_executable), str(mtz_path)],
            env=dict(environment) if environment is not None else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ReflectionDataError(f"Could not inspect MTZ columns: {exc}") from exc
    if completed.returncode:
        raise ReflectionDataError(
            f"phenix.mtz.dump failed ({completed.returncode}) while listing "
            f"the columns of {mtz_path.name}"
        )
    return completed.stdout


def resolve_phaser_data(
    choice: str,
    mtz_path: Path,
    mtz_dump_executable: Path | None,
    environment: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Freeze which reflection array Phaser will read.

    ``auto`` leaves ``phaser.labin`` unset, so Phenix keeps its own choice
    (amplitudes when present). ``intensities`` selects the labels here so the
    preflight fails before Phaser if the MTZ has no merged intensities.
    """
    selected = normalize_phaser_data(choice)
    if selected == AUTO:
        return {
            "schema_version": SCHEMA_VERSION,
            "requested": AUTO,
            "labin": None,
            "selection": "phenix.phaser default (amplitudes when present)",
        }
    if mtz_dump_executable is None:
        raise ReflectionDataError(
            "phaser_data = intensities requires phenix.mtz.dump to read the MTZ "
            "column types"
        )
    columns = parse_mtz_dump_columns(
        read_mtz_dump(mtz_path, mtz_dump_executable, environment)
    )
    labin = select_intensity_labels(columns)
    return {
        "schema_version": SCHEMA_VERSION,
        "requested": INTENSITIES,
        "labin": labin,
        "selection": "first merged intensity array (MTZ column types J,Q)",
        "mtz_columns": [f"{column.label}:{column.type}" for column in columns],
    }


_LABIN_RE = re.compile(r"^\s*LABIN\s+(.+?)\s*$", re.M)


def reported_labin(phaser_output: str) -> str | None:
    """Return Phaser's own ``LABIN`` echo, for example ``I = I SIGI = SIGI``."""
    match = _LABIN_RE.search(phaser_output)
    return match.group(1) if match else None


def reported_data_kind(labin_text: str | None) -> str | None:
    """Classify a Phaser ``LABIN`` echo as intensities or amplitudes."""
    if not labin_text:
        return None
    first = labin_text.split("=", 1)[0].strip().upper()
    if first in {"I", "IMEAN", "I(+)"}:
        return INTENSITIES
    if first in {"F", "FMEAN", "F(+)"}:
        return "amplitudes"
    return None


__all__ = [
    "AUTO",
    "INTENSITIES",
    "MtzColumn",
    "PHASER_DATA_CHOICES",
    "ReflectionDataError",
    "mean_intensity_label_strings",
    "normalize_phaser_data",
    "parse_mtz_dump_columns",
    "reported_data_kind",
    "reported_labin",
    "resolve_phaser_data",
    "select_intensity_labels",
]
