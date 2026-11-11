"""Portable paths for committed NASolve datasets and run artifacts.

Operational references are anchored to the run, dataset, or repository instead
of a collaborator's absolute checkout path.  The resolver also understands the
absolute paths written by schema-1 run reports so historical runs remain usable
after a clone or checkout move.
"""

from __future__ import annotations

import hashlib
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Mapping


RUN_ANCHOR = "run"
DATASET_ANCHOR = "dataset"
REPOSITORY_ANCHOR = "repository"
ABSOLUTE_ANCHOR = "absolute"


# Numbered run series that live directly in a dataset directory.
_DATASET_RUN_SERIES = {"AutoMR", "AutoSolMetal"}


def dataset_directory(run_directory: Path) -> Path:
    run = run_directory.expanduser().resolve()
    if run.parent.name in _DATASET_RUN_SERIES:
        return run.parent.parent
    return run.parent


def repository_directory(run_directory: Path) -> Path | None:
    run = run_directory.expanduser().resolve()
    for candidate in (run, *run.parents):
        if (
            (candidate / "pyproject.toml").is_file()
            and (candidate / "MR_frames").is_dir()
        ):
            return candidate
    return None


def artifact_reference(path: Path, run_directory: Path) -> dict[str, str]:
    """Return a schema-2 portable reference for an existing artifact."""
    run = run_directory.expanduser().resolve()
    dataset = dataset_directory(run)
    repository = repository_directory(run)
    resolved = path.expanduser().resolve()
    anchors = [
        (RUN_ANCHOR, run),
        (DATASET_ANCHOR, dataset),
    ]
    if repository is not None:
        anchors.append((REPOSITORY_ANCHOR, repository))
    for anchor, root in anchors:
        try:
            relative = resolved.relative_to(root)
        except ValueError:
            continue
        return {"anchor": anchor, "relative_path": relative.as_posix()}
    return {"anchor": ABSOLUTE_ANCHOR, "absolute_path": str(resolved)}


def _anchored_candidate(path: Path, anchor: str, run: Path) -> Path | None:
    roots = {
        RUN_ANCHOR: run,
        DATASET_ANCHOR: dataset_directory(run),
        REPOSITORY_ANCHOR: repository_directory(run),
    }
    if anchor == ABSOLUTE_ANCHOR:
        if not path.is_absolute() or ".." in path.parts:
            return None
        return path.expanduser().resolve()
    root = roots.get(anchor)
    if root is None or path.is_absolute():
        return None
    candidate = (root / path).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError:
        return None
    return candidate


def _safe_join(root: Path, parts: tuple[str, ...]) -> Path | None:
    if not parts or any(part in {"", ".", ".."} for part in parts):
        return None
    resolved_root = root.resolve()
    candidate = resolved_root.joinpath(*parts).resolve()
    try:
        candidate.relative_to(resolved_root)
    except ValueError:
        return None
    return candidate


def _legacy_parts(value: str) -> list[tuple[str, ...]]:
    variants = [PurePosixPath(value).parts, PureWindowsPath(value).parts]
    unique: list[tuple[str, ...]] = []
    for parts in variants:
        if parts not in unique:
            unique.append(parts)
    return unique


def _legacy_rebased_candidates(value: str, run: Path) -> list[Path]:
    dataset = dataset_directory(run)
    repository = repository_directory(run)
    candidates: list[Path] = []
    path = Path(value).expanduser()
    windows_path = PureWindowsPath(value)
    if not path.is_absolute() and not windows_path.is_absolute():
        for root in (run, dataset):
            candidate = _safe_join(root, path.parts)
            if candidate is not None:
                candidates.append(candidate)

    for parts in _legacy_parts(value):
        folded = [part.casefold() for part in parts]
        for index in range(len(parts) - 1):
            if (
                folded[index] == "automr"
                and folded[index + 1] == run.name.casefold()
            ):
                candidate = _safe_join(run, parts[index + 2 :])
                if candidate is not None:
                    candidates.append(candidate)
        for index, part in enumerate(folded):
            if part == dataset.name.casefold():
                candidate = _safe_join(dataset, parts[index + 1 :])
                if candidate is not None:
                    candidates.append(candidate)
            if part == "mr_frames" and repository is not None:
                candidate = _safe_join(repository, parts[index:])
                if candidate is not None:
                    candidates.append(candidate)
    return candidates


def _checksum_matches(path: Path, expected: object) -> bool:
    if not isinstance(expected, str) or not expected:
        return True
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest() == expected.casefold()


def resolve_artifact_path(
    value: object,
    run_directory: Path,
    *,
    anchor: str | None = None,
    verify_checksum: bool = True,
) -> Path | None:
    """Resolve a portable or legacy artifact path if its local file exists."""
    expected_checksum: object = None
    if isinstance(value, Mapping):
        anchor_value = value.get("anchor")
        anchor = anchor_value if isinstance(anchor_value, str) else anchor
        expected_checksum = value.get("sha256")
        if anchor == ABSOLUTE_ANCHOR:
            value = value.get("absolute_path", value.get("path"))
        else:
            value = value.get("relative_path", value.get("path"))
    if not isinstance(value, (str, Path)) or not str(value):
        return None
    run = run_directory.expanduser().resolve()
    text_value = str(value)
    path = Path(value).expanduser()
    if anchor is not None:
        candidate = _anchored_candidate(path, anchor, run)
        if candidate is None or not candidate.is_file():
            return None
        if verify_checksum and not _checksum_matches(candidate, expected_checksum):
            return None
        return candidate

    parts_variants = _legacy_parts(text_value)
    if any(".." in parts for parts in parts_variants):
        return None
    matches: list[Path] = []
    local_candidate_exists = False
    for candidate in _legacy_rebased_candidates(text_value, run):
        resolved = candidate.resolve()
        if not resolved.is_file():
            continue
        local_candidate_exists = True
        if resolved in matches:
            continue
        if verify_checksum and not _checksum_matches(resolved, expected_checksum):
            continue
        matches.append(resolved)
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        return None
    if local_candidate_exists:
        return None

    # An unchanged legacy absolute path is the final compatibility fallback.
    if path.is_absolute() and path.is_file():
        resolved = path.resolve()
        if not verify_checksum or _checksum_matches(resolved, expected_checksum):
            return resolved
    return None


__all__ = [
    "ABSOLUTE_ANCHOR",
    "DATASET_ANCHOR",
    "REPOSITORY_ANCHOR",
    "RUN_ANCHOR",
    "artifact_reference",
    "dataset_directory",
    "repository_directory",
    "resolve_artifact_path",
]
