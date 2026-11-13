"""Choose which observed reflection array Phaser reads from the input MTZ.

AutoMR always tells Phaser which data to use: ``F`` (merged amplitudes, the
default) or ``I`` (merged intensities, what the Phenix GUI uses and what
Phaser's intensity likelihood target is designed for). Nothing is left to
``phenix.phaser``'s own fallback, so the choice is visible in ``nasolve.txt``,
the run snapshot, the report and ``phaser.eff``.

``phaser.labin`` must equal the label string that ``iotbx`` gives the array,
for example ``F,SIGF`` or ``IMEAN,SIGIMEAN``. The labels are therefore chosen
from the MTZ column *types* reported by ``phenix.mtz.dump`` (``F``/``J``
followed by ``Q``), grouped the way ``iotbx.mtz.object.group_columns`` groups
them, rather than from guessed column names.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping


AMPLITUDES = "F"
INTENSITIES = "I"
DEFAULT_PHASER_DATA = AMPLITUDES
PHASER_DATA_CHOICES = (AMPLITUDES, INTENSITIES)
_ALIASES = {
    "f": AMPLITUDES,
    "amplitude": AMPLITUDES,
    "amplitudes": AMPLITUDES,
    "fsigf": AMPLITUDES,
    "f/sigf": AMPLITUDES,
    "f,sigf": AMPLITUDES,
    "i": INTENSITIES,
    "intensity": INTENSITIES,
    "intensities": INTENSITIES,
    "isigi": INTENSITIES,
    "i/sigi": INTENSITIES,
    "i,sigi": INTENSITIES,
}
_DESCRIPTION = {AMPLITUDES: "amplitudes", INTENSITIES: "intensities"}
SCHEMA_VERSION = 2
# phenix.phaser's own preference among amplitude arrays (driver.run).
_PREFERRED_AMPLITUDES = ("F,SIGF", "FOBS,SIGFOBS")


class ReflectionDataError(RuntimeError):
    """Raised when the requested Phaser data cannot be selected safely."""


@dataclass(frozen=True)
class MtzColumn:
    dataset: int  # ordinal of the MTZ dataset in phenix.mtz.dump order
    label: str
    type: str
    valid: int | None = None
    minimum: float | None = None
    maximum: float | None = None


@dataclass(frozen=True)
class MtzDataset:
    ordinal: int  # position in phenix.mtz.dump order
    name: str | None
    wavelength: float | None
    columns: tuple[MtzColumn, ...]


def describe(choice: str) -> str:
    return _DESCRIPTION.get(choice, choice)


def normalize_phaser_data(value: object) -> str:
    """Return ``F`` or ``I`` for a user setting (default ``F``)."""
    if value is None:
        return DEFAULT_PHASER_DATA
    if not isinstance(value, str):
        raise ReflectionDataError("phaser_data must be F or I")
    key = "".join(value.split()).casefold()
    if key not in _ALIASES:
        raise ReflectionDataError(
            f"phaser_data must be F (amplitudes) or I (intensities), not {value!r}"
        )
    return _ALIASES[key]


# A column row of ``phenix.mtz.dump``'s summary, for example
# ``    IMEAN   1635 100.00%  -93.85  4926.23 J: intensity``
_COLUMN_RE = re.compile(
    r"^\s+(\S+)\s+(\d+)\s+[0-9.]+%\s+(\S+)\s+(\S+)\s+([A-Za-z]):(?:\s|$)"
)
_DATASET_RE = re.compile(r"^\s*Dataset\s+(\d+)\s*:\s*$")
_CRYSTAL_RE = re.compile(r"^\s*Crystal\s+(\d+)\s*:\s*$")
_DATASET_NAME_RE = re.compile(r"^\s+Name:\s*(.*?)\s*$")
_WAVELENGTH_RE = re.compile(r"^\s+Wavelength:\s*(\S+)\s*$")


def _number(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


def parse_mtz_dump_columns(text: str) -> tuple[MtzColumn, ...]:
    """Return the columns listed by ``phenix.mtz.dump``, in file order."""
    columns: list[MtzColumn] = []
    # Dataset numbers restart in every crystal, so number them in file order.
    dataset = 0
    for line in text.splitlines():
        if _DATASET_RE.match(line):
            dataset += 1
            continue
        match = _COLUMN_RE.match(line)
        if match:
            label, valid, minimum, maximum, kind = match.groups()
            columns.append(MtzColumn(
                dataset, label, kind.upper(), int(valid), _number(minimum), _number(maximum)
            ))
    return tuple(columns)


def parse_mtz_dump_datasets(text: str) -> tuple[MtzDataset, ...]:
    """Return each MTZ dataset's name, wavelength and columns."""
    datasets: list[MtzDataset] = []
    current: dict[str, object] | None = None

    def close() -> None:
        if current is not None:
            datasets.append(MtzDataset(
                len(datasets) + 1, current["name"], current["wavelength"],
                tuple(current["columns"]),
            ))

    ordinal = 0
    for line in text.splitlines():
        if _CRYSTAL_RE.match(line):
            close()
            current = None
            continue
        if _DATASET_RE.match(line):
            close()
            ordinal += 1
            current = {"name": None, "wavelength": None, "columns": []}
            continue
        if current is None:
            continue
        name = _DATASET_NAME_RE.match(line)
        if name and current["name"] is None and not current["columns"]:
            current["name"] = name.group(1)
            continue
        wavelength = _WAVELENGTH_RE.match(line)
        if wavelength:
            current["wavelength"] = _number(wavelength.group(1))
            continue
        match = _COLUMN_RE.match(line)
        if match:
            label, valid, minimum, maximum, kind = match.groups()
            current["columns"].append(MtzColumn(
                ordinal, label, kind.upper(), int(valid), _number(minimum), _number(maximum)
            ))
    close()
    return tuple(datasets)


def mtz_wavelength(text: str) -> tuple[float, str | None]:
    """Return the X-ray wavelength recorded in the MTZ header.

    Only datasets that hold reflection data count: the ``HKL_base`` dataset
    carries just H, K, L and a placeholder wavelength. Every data dataset must
    report the same plausible value.
    """
    data_sets = [
        dataset for dataset in parse_mtz_dump_datasets(text)
        if any(column.type != "H" for column in dataset.columns)
    ]
    if not data_sets:
        raise ReflectionDataError(
            "phenix.mtz.dump listed no MTZ dataset with reflection data, so the "
            "wavelength could not be read from the MTZ header"
        )
    for dataset in data_sets:
        if dataset.wavelength is None:
            raise ReflectionDataError(
                f"MTZ dataset {dataset.name!r} has no wavelength in its header"
            )
        if not 0.1 <= dataset.wavelength <= 5.0:
            raise ReflectionDataError(
                f"MTZ dataset {dataset.name!r} has an implausible wavelength "
                f"{dataset.wavelength:g} A in its header"
            )
    reference = data_sets[0].wavelength
    assert reference is not None
    if any(abs(float(dataset.wavelength) - reference) > 5e-5 for dataset in data_sets):
        values = ", ".join(f"{dataset.name}={dataset.wavelength:g}" for dataset in data_sets)
        raise ReflectionDataError(f"MTZ datasets disagree on the wavelength: {values}")
    return reference, data_sets[0].name


def _is_anomalous_label(label: str) -> bool:
    upper = label.upper()
    return (
        "(+)" in upper or "(-)" in upper
        or upper.endswith(("+", "-", "PLUS", "MINUS"))
    )


def _isym_label(columns: list[MtzColumn], start: int) -> str | None:
    """The ISYM column iotbx appends to an F,SIGF,DANO,SIGDANO group."""
    for column in columns[start:]:
        if column.type == "Y":
            if (
                column.valid and column.minimum is not None and column.maximum is not None
                and column.minimum >= 0 and column.maximum <= 2
            ):
                return column.label
            return None
    return None


def merged_arrays(columns: Iterable[MtzColumn]) -> list[tuple[str, str]]:
    """Return ``(F|I, iotbx label string)`` for merged observation arrays.

    Walks each MTZ dataset in file order like ``iotbx``: ``F Q D Q`` is one
    reconstructed-amplitude array (plus the following ISYM), ``F``/``J``
    followed by ``Q`` is one array, and anomalous ``G L``/``K M`` groups are
    skipped. A column is therefore only paired with the column ``iotbx``
    itself pairs it with.
    """
    by_dataset: dict[int, list[MtzColumn]] = {}
    for column in columns:
        by_dataset.setdefault(column.dataset, []).append(column)
    found: list[tuple[str, str]] = []
    for dataset_columns in by_dataset.values():
        types = "".join(column.type for column in dataset_columns)
        index = 0
        while index < len(dataset_columns):
            kind = types[index]
            if types[index:index + 4] == "FQDQ":
                labels = [column.label for column in dataset_columns[index:index + 4]]
                isym = _isym_label(dataset_columns, index + 4)
                if isym is not None:
                    labels.append(isym)
                found.append((AMPLITUDES, ",".join(labels)))
                index += 4
                continue
            if kind in "JFED" and index + 1 < len(types) and types[index + 1] in "QP":
                first, second = dataset_columns[index], dataset_columns[index + 1]
                if (
                    kind in "JF" and types[index + 1] == "Q"
                    and not _is_anomalous_label(first.label)
                    and not _is_anomalous_label(second.label)
                ):
                    found.append((
                        INTENSITIES if kind == "J" else AMPLITUDES,
                        f"{first.label},{second.label}",
                    ))
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


def select_labels(columns: Iterable[MtzColumn], choice: str) -> str:
    """Return the ``labin`` for ``F`` or ``I``, or explain why there is none.

    ``I``: the first merged intensity array. ``F``: ``F,SIGF`` or
    ``FOBS,SIGFOBS`` when present, otherwise the first merged amplitude array
    that is not a calculated ``FC`` column (phenix.phaser's own preference).
    """
    columns = tuple(columns)
    if not columns:
        raise ReflectionDataError(
            "phenix.mtz.dump did not list any MTZ columns, so the Phaser data "
            "columns could not be identified"
        )
    candidates = [labels for kind, labels in merged_arrays(columns) if kind == choice]
    if choice == AMPLITUDES:
        preferred = [labels for labels in candidates if labels in _PREFERRED_AMPLITUDES]
        candidates = preferred or [
            labels for labels in candidates if not labels.upper().startswith("FC")
        ]
    if not candidates:
        present = ", ".join(f"{column.label}({column.type})" for column in columns)
        needed = (
            "a merged amplitude column (MTZ type F)" if choice == AMPLITUDES
            else "a merged intensity column (MTZ type J)"
        )
        raise ReflectionDataError(
            f"phaser_data = {choice} needs {needed} followed by its standard "
            "deviation (type Q); the MTZ has: " + present
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
    """Freeze which reflection array Phaser will read (always explicit)."""
    selected = normalize_phaser_data(choice)
    if mtz_dump_executable is None:
        raise ReflectionDataError(
            "Choosing the Phaser data columns requires phenix.mtz.dump to read "
            "the MTZ column types"
        )
    columns = parse_mtz_dump_columns(
        read_mtz_dump(mtz_path, mtz_dump_executable, environment)
    )
    labin = select_labels(columns, selected)
    return {
        "schema_version": SCHEMA_VERSION,
        "requested": selected,
        "labin": labin,
        "selection": (
            "first merged intensity array (MTZ column types J,Q)"
            if selected == INTENSITIES else
            "F,SIGF or FOBS,SIGFOBS, else first merged non-FC amplitude array"
        ),
        "mtz_columns": [f"{column.label}:{column.type}" for column in columns],
    }


_LABIN_RE = re.compile(r"^\s*LABIN\s+(.+?)\s*$", re.M)


def reported_labin(phaser_output: str) -> str | None:
    """Return Phaser's own ``LABIN`` echo, for example ``I = I SIGI = SIGI``."""
    match = _LABIN_RE.search(phaser_output)
    return match.group(1) if match else None


def reported_data_kind(labin_text: str | None) -> str | None:
    """Classify a Phaser ``LABIN`` echo as ``I`` or ``F``."""
    if not labin_text:
        return None
    first = labin_text.split("=", 1)[0].strip().upper()
    if first in {"I", "IMEAN", "I(+)"}:
        return INTENSITIES
    if first in {"F", "FMEAN", "F(+)"}:
        return AMPLITUDES
    return None


__all__ = [
    "AMPLITUDES",
    "DEFAULT_PHASER_DATA",
    "INTENSITIES",
    "MtzColumn",
    "MtzDataset",
    "PHASER_DATA_CHOICES",
    "ReflectionDataError",
    "describe",
    "merged_arrays",
    "mtz_wavelength",
    "parse_mtz_dump_datasets",
    "normalize_phaser_data",
    "parse_mtz_dump_columns",
    "reported_data_kind",
    "reported_labin",
    "resolve_phaser_data",
    "select_labels",
]
