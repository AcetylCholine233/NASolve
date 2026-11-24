"""Immutable, user-visible refinement checkpoints.

The checkpoint registry is deliberately independent of Phenix.  It records the
model lineage and the frozen crystallographic inputs needed to reproduce or
branch a refinement, including work performed manually outside NASolve.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

from .phosphate import PhosphateError
from .ligand_profiles import validate_model_phosphate_policy
from .run_context import artifact_reference, resolve_artifact_path


class CheckpointError(RuntimeError):
    """Raised when checkpoint provenance cannot be established safely."""


@dataclass(frozen=True)
class CheckpointRecord:
    checkpoint_id: str
    parent: str | None
    status: str
    model: Path | None
    usable: bool
    recipe: str
    label: str | None
    metrics: dict[str, float | None]


SCHEMA_VERSION = 2
READABLE_SCHEMA_VERSIONS = {1, SCHEMA_VERSION}
AUTOMATIC_STATUSES = {"READY", "SUCCESS", "USER_APPROVED"}
REUSABLE_STATUSES = AUTOMATIC_STATUSES | {"REVIEW"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _required_file(value: object, description: str, run: Path) -> Path:
    if not isinstance(value, str) or not value:
        raise CheckpointError(f"Run report has no {description}")
    path = resolve_artifact_path(value, run)
    if path is None:
        raise CheckpointError(f"Required {description} is missing: {value}")
    return path


def _run_and_report(run_directory: Path) -> tuple[Path, dict[str, object]]:
    selected = run_directory.expanduser().resolve()
    report_path = selected if selected.name == "report.json" else selected / "report.json"
    run = report_path.parent
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CheckpointError(f"Could not read NASolve report {report_path}: {exc}") from exc
    if not isinstance(report, dict):
        raise CheckpointError(f"NASolve report is not a JSON object: {report_path}")
    if report.get("workflow") != "automr":
        raise CheckpointError("Checkpoints require a NASolve AutoMR run")
    postmr = report.get("postmr")
    if not isinstance(postmr, Mapping) or postmr.get("status") != "POSTMR_READY":
        raise CheckpointError("Checkpoints require a completed POSTMR_READY model")
    return run, report


def registry_path(run_directory: Path) -> Path:
    run = run_directory.expanduser().resolve()
    if run.name == "report.json":
        run = run.parent
    return run / "AutoRefine" / "checkpoints.json"


def _file_reference(
    path: Path, run: Path, *, sha256: str | None = None
) -> dict[str, object]:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise CheckpointError(f"Checkpoint input is missing: {resolved}")
    return {
        **artifact_reference(resolved, run),
        "sha256": sha256 or _sha256(resolved),
        "size": resolved.stat().st_size,
    }


def _is_portable_reference(value: Mapping[str, object]) -> bool:
    return isinstance(value.get("anchor"), str) and (
        isinstance(value.get("relative_path"), str)
        or isinstance(value.get("absolute_path"), str)
    )


def _require_portable_checksum(value: Mapping[str, object], description: str) -> None:
    expected = value.get("sha256")
    if (
        not isinstance(expected, str)
        or re.fullmatch(r"[0-9a-fA-F]{64}", expected) is None
    ):
        raise CheckpointError(f"Checkpoint {description} has a malformed checksum")


def _legacy_provenance(value: object) -> dict[str, object] | None:
    if isinstance(value, str):
        return {"provenance_path": value, "available": False}
    if not isinstance(value, Mapping):
        return None
    migrated = dict(value)
    legacy_path = migrated.pop("path", None)
    if isinstance(legacy_path, str):
        migrated.setdefault("provenance_path", legacy_path)
    migrated["available"] = False
    return migrated


def _ref_path(
    value: object,
    description: str,
    run: Path | None,
    *,
    required: bool = True,
    verify_checksum: bool = True,
) -> Path | None:
    path_value = (
        value.get("relative_path", value.get("absolute_path", value.get("path")))
        if isinstance(value, Mapping)
        else None
    )
    if not isinstance(path_value, str):
        if required:
            raise CheckpointError(f"Checkpoint has no {description}")
        return None
    portable = isinstance(value, Mapping) and _is_portable_reference(value)
    expected = value.get("sha256") if isinstance(value, Mapping) else None
    if portable:
        assert isinstance(value, Mapping)
        _require_portable_checksum(value, description)
    if run is not None:
        path = resolve_artifact_path(value, run, verify_checksum=False)
    elif isinstance(value, Mapping) and isinstance(value.get("anchor"), str):
        path = None
    else:
        path = Path(path_value).expanduser().resolve()
    if path is None or not path.is_file():
        if required:
            raise CheckpointError(f"Checkpoint {description} is missing: {path_value}")
        return None
    if (
        verify_checksum
        and isinstance(expected, str)
        and _sha256(path) != expected.casefold()
    ):
        raise CheckpointError(
            f"Checkpoint {description} changed after it was frozen: {path}"
        )
    return path


def _initial_restraints(postmr: Mapping[str, object], run: Path) -> list[Path]:
    if "refinement_restraints" in postmr:
        references = postmr["refinement_restraints"]
        if not isinstance(references, list) or not references:
            raise CheckpointError("Malformed frozen refinement restraint list")
        return [_ref_path(ref, "PostMR restraint", run) for ref in references]
    values = postmr.get("restraints")
    reported = (
        [value for value in values if isinstance(value, str)]
        if isinstance(values, list)
        else []
    )
    non_cif_values = [
        value for value in reported if not value.casefold().endswith(".cif")
    ]
    cif_values = [value for value in reported if value.casefold().endswith(".cif")]
    retained: list[Path] = []
    for value in non_cif_values:
        resolved = resolve_artifact_path(value, run)
        if resolved is None:
            raise CheckpointError(f"PostMR restraint input is missing: {value}")
        retained.append(resolved)
    readyset = postmr.get("readyset")
    generated = readyset.get("generated_ligand_cif") if isinstance(readyset, Mapping) else None
    generated_path = resolve_artifact_path(generated, run)
    # ReadySet's combined ligand dictionary supersedes the individual input CIFs.
    if generated_path is not None:
        retained.append(generated_path)
    else:
        for value in cif_values:
            resolved = resolve_artifact_path(value, run)
            if resolved is None:
                raise CheckpointError(f"PostMR restraint input is missing: {value}")
            retained.append(resolved)
    return list(dict.fromkeys(retained))


def _phase_reference(report: Mapping[str, object], run: Path) -> dict[str, object] | None:
    autosol = report.get("autosol")
    if not isinstance(autosol, Mapping) or autosol.get("use_for_refinement") is not True:
        return None
    outputs = autosol.get("outputs")
    if not isinstance(outputs, Mapping):
        raise CheckpointError("Approved AutoSol result has no output references")
    path = _required_file(outputs.get("refinement_data"), "AutoSol refinement data", run)
    actual_sha256 = _sha256(path)
    expected = outputs.get("refinement_data_sha256")
    if expected is not None and (
        not isinstance(expected, str)
        or re.fullmatch(r"[0-9a-fA-F]{64}", expected) is None
    ):
        raise CheckpointError("Approved AutoSol refinement-data checksum is malformed")
    if isinstance(expected, str) and actual_sha256 != expected.casefold():
        raise CheckpointError(
            "Approved AutoSol refinement data failed checksum validation"
        )
    labels = outputs.get("refinement_data_phase_labels")
    return {
        **_file_reference(path, run, sha256=actual_sha256),
        "labels": list(labels) if isinstance(labels, list) else [],
        "source": "validated-autosol",
    }


def _root_payload(report: Mapping[str, object], run: Path) -> dict[str, object]:
    postmr = report["postmr"]
    assert isinstance(postmr, Mapping)
    inputs = report.get("inputs")
    if not isinstance(inputs, Mapping):
        raise CheckpointError("Run report has no frozen AutoMR inputs")
    model = _required_file(postmr.get("prepared_model"), "PostMR prepared model", run)
    try:
        validate_model_phosphate_policy(model, report)
    except PhosphateError as exc:
        raise CheckpointError(str(exc)) from exc
    model_sha256 = _sha256(model)
    expected_model_sha256 = postmr.get("prepared_sha256")
    if expected_model_sha256 is not None and not isinstance(
        expected_model_sha256, str
    ):
        raise CheckpointError("PostMR prepared-model checksum is malformed")
    model_verified = isinstance(expected_model_sha256, str)
    if (
        isinstance(expected_model_sha256, str)
        and model_sha256 != expected_model_sha256.casefold()
    ):
        raise CheckpointError(
            "Rebased PostMR prepared model does not match its frozen checksum"
        )
    reflections = _required_file(
        inputs.get("reflections"), "authoritative reflections", run
    )
    reflections_sha256 = _sha256(reflections)
    expected_reflections_sha256 = inputs.get("reflections_sha256")
    if expected_reflections_sha256 is not None and not isinstance(
        expected_reflections_sha256, str
    ):
        raise CheckpointError("Authoritative reflections checksum is malformed")
    observations_verified = isinstance(expected_reflections_sha256, str)
    if (
        observations_verified
        and reflections_sha256 != expected_reflections_sha256.casefold()
    ):
        raise CheckpointError(
            "Rebased authoritative reflections do not match their frozen checksum"
        )
    restraint_refs = [
        _file_reference(path, run) for path in _initial_restraints(postmr, run)
    ]
    return {
        "id": "postmr",
        "parent": None,
        "kind": "postmr",
        "status": "READY",
        "usable": True,
        "created_utc": str(postmr.get("created_utc") or _now()),
        "recipe": "PostMR/ReadySet",
        "label": "ReadySet model",
        "model": _file_reference(model, run, sha256=model_sha256),
        "observations": _file_reference(
            reflections, run, sha256=reflections_sha256
        ),
        "phases": _phase_reference(report, run),
        "restraints": restraint_refs,
        "metrics": {},
        "compatibility": {
            "validated": True,
            "source": "POSTMR_READY",
            "model_checksum": (
                "verified" if model_verified else "legacy-unverified"
            ),
            "observations_checksum": (
                "verified" if observations_verified else "legacy-unverified"
            ),
        },
    }


def read_registry(run_directory: Path) -> tuple[Path, dict[str, object]]:
    """Read a registry or derive its validated PostMR root without writing files."""
    run, report = _run_and_report(run_directory)
    path = registry_path(run)
    if path.is_file():
        try:
            registry = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise CheckpointError(f"Could not read checkpoint registry {path}: {exc}") from exc
        if not isinstance(registry, dict):
            raise CheckpointError(f"Checkpoint registry is not a JSON object: {path}")
        schema_version = registry.get("schema_version")
        if (
            type(schema_version) is not int
            or schema_version not in READABLE_SCHEMA_VERSIONS
        ):
            raise CheckpointError(f"Unsupported checkpoint schema in {path}")
        checkpoints = registry.get("checkpoints")
        if not isinstance(checkpoints, list) or not any(
            isinstance(item, Mapping) and item.get("id") == "postmr" for item in checkpoints
        ):
            raise CheckpointError("Checkpoint registry has no PostMR root")
        return run, registry
    registry = {
        "schema_version": SCHEMA_VERSION,
        "run_directory": ".",
        "created_utc": _now(),
        "updated_utc": _now(),
        "current": "postmr",
        "bookmarks": {},
        "checkpoints": [_root_payload(report, run)],
    }
    return run, registry


def initialize_registry(run_directory: Path) -> tuple[Path, dict[str, object]]:
    """Persist the root PostMR checkpoint when a modifying operation needs it."""
    run, registry = read_registry(run_directory)
    path = registry_path(run)
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        _write_json(path, registry)
    return run, registry


def save_registry(run: Path, registry: dict[str, object]) -> Path:
    registry["schema_version"] = SCHEMA_VERSION
    registry["run_directory"] = "."
    for item in _items(registry):
        for key in ("model", "observations", "phases"):
            value = item.get(key)
            if not isinstance(value, dict):
                continue
            if _is_portable_reference(value):
                _require_portable_checksum(value, key)
                continue
            path = _ref_path(value, key, run, required=False)
            if path is not None:
                portable = artifact_reference(path, run)
                if portable.get("anchor") != "absolute":
                    legacy_path = value.pop("path", None)
                    if isinstance(legacy_path, str):
                        value.setdefault("provenance_path", legacy_path)
                    value.update(portable)
                    continue
            provenance = _legacy_provenance(value)
            if provenance is not None:
                item[key] = provenance
        restraints = item.get("restraints")
        if isinstance(restraints, list):
            for value in restraints:
                if not isinstance(value, dict):
                    continue
                if _is_portable_reference(value):
                    _require_portable_checksum(value, "restraint")
                    continue
                path = _ref_path(value, "restraint", run, required=False)
                if path is not None:
                    portable = artifact_reference(path, run)
                    if portable.get("anchor") != "absolute":
                        legacy_path = value.pop("path", None)
                        if isinstance(legacy_path, str):
                            value.setdefault("provenance_path", legacy_path)
                        value.update(portable)
                        continue
                provenance = _legacy_provenance(value)
                if provenance is not None:
                    value.clear()
                    value.update(provenance)
        outputs = item.get("outputs")
        if isinstance(outputs, dict):
            for key, value in tuple(outputs.items()):
                if isinstance(value, Mapping) and _is_portable_reference(value):
                    _require_portable_checksum(value, f"output {key}")
                    continue
                path = resolve_artifact_path(value, run)
                if path is not None:
                    reference = _file_reference(path, run)
                    if reference.get("anchor") != "absolute":
                        outputs[key] = reference
                        continue
                provenance = _legacy_provenance(value)
                if provenance is not None:
                    outputs[key] = provenance
    registry["updated_utc"] = _now()
    path = registry_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_json(path, registry)
    return path


def _items(registry: Mapping[str, object]) -> list[dict[str, object]]:
    values = registry.get("checkpoints")
    if not isinstance(values, list) or not all(isinstance(item, dict) for item in values):
        raise CheckpointError("Malformed checkpoint registry")
    return values  # type: ignore[return-value]


def resolve_checkpoint(
    registry: Mapping[str, object], checkpoint: str | None = None
) -> dict[str, object]:
    requested = checkpoint or registry.get("current")
    if not isinstance(requested, str) or not requested:
        raise CheckpointError("No current checkpoint is selected")
    bookmarks = registry.get("bookmarks")
    if isinstance(bookmarks, Mapping) and isinstance(bookmarks.get(requested), str):
        requested = bookmarks[requested]
    for item in _items(registry):
        if item.get("id") == requested:
            return item
    raise CheckpointError(f"Unknown checkpoint or bookmark: {requested}")


def checkpoint_record(
    item: Mapping[str, object], run: Path | None = None
) -> CheckpointRecord:
    metrics = item.get("metrics")
    normalized_metrics = {
        str(key): float(value) if isinstance(value, (int, float)) else None
        for key, value in metrics.items()
    } if isinstance(metrics, Mapping) else {}
    return CheckpointRecord(
        checkpoint_id=str(item.get("id")),
        parent=str(item["parent"]) if isinstance(item.get("parent"), str) else None,
        status=str(item.get("status", "UNKNOWN")),
        model=_ref_path(item.get("model"), "model", run, required=False),
        usable=item.get("usable") is True,
        recipe=str(item.get("recipe", "unknown")),
        label=str(item["label"]) if isinstance(item.get("label"), str) else None,
        metrics=normalized_metrics,
    )


def list_checkpoints(run_directory: Path) -> tuple[list[CheckpointRecord], str, dict[str, str]]:
    run, registry = read_registry(run_directory)
    current = str(registry.get("current", "postmr"))
    bookmarks = registry.get("bookmarks")
    aliases = {
        str(name): str(target)
        for name, target in bookmarks.items()
        if isinstance(bookmarks, Mapping) and isinstance(name, str) and isinstance(target, str)
    } if isinstance(bookmarks, Mapping) else {}
    return [checkpoint_record(item, run) for item in _items(registry)], current, aliases


def next_checkpoint_id(registry: Mapping[str, object], prefix: str) -> str:
    expression = re.compile(rf"^{re.escape(prefix)}-(\d+)$")
    numbers = [
        int(match.group(1))
        for item in _items(registry)
        if isinstance(item.get("id"), str)
        and (match := expression.match(str(item["id"]))) is not None
    ]
    return f"{prefix}-{max(numbers, default=0) + 1:03d}"


def append_checkpoint(
    run: Path,
    registry: dict[str, object],
    payload: dict[str, object],
    *,
    select: bool = False,
) -> Path:
    checkpoint_id = payload.get("id")
    if not isinstance(checkpoint_id, str) or not checkpoint_id:
        raise CheckpointError("New checkpoint has no ID")
    if any(item.get("id") == checkpoint_id for item in _items(registry)):
        raise CheckpointError(f"Checkpoint already exists: {checkpoint_id}")
    payload.setdefault("created_utc", _now())
    _items(registry).append(payload)
    if select:
        registry["current"] = checkpoint_id
    return save_registry(run, registry)


def select_checkpoint(run_directory: Path, checkpoint: str) -> CheckpointRecord:
    run, registry = initialize_registry(run_directory)
    item = resolve_checkpoint(registry, checkpoint)
    if item.get("usable") is not True or item.get("status") not in REUSABLE_STATUSES:
        raise CheckpointError(f"Checkpoint {item.get('id')} is not reusable")
    _ref_path(item.get("model"), "model", run)
    registry["current"] = item["id"]
    save_registry(run, registry)
    return checkpoint_record(item, run)


def _valid_bookmark(name: str) -> str:
    normalized = name.strip()
    if not normalized or len(normalized) > 80:
        raise CheckpointError("Checkpoint name must contain 1-80 characters")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_. -]*", normalized):
        raise CheckpointError(
            "Checkpoint name may contain letters, numbers, spaces, '.', '_' and '-'"
        )
    return normalized


def add_checkpoint(
    run_directory: Path,
    *,
    name: str,
    model: Path | None = None,
    reflections: Path | None = None,
    parent: str | None = None,
) -> CheckpointRecord:
    """Bookmark the current node, or import a manual model as a new child."""
    run, registry = initialize_registry(run_directory)
    label = _valid_bookmark(name)
    bookmarks = registry.setdefault("bookmarks", {})
    if not isinstance(bookmarks, dict):
        raise CheckpointError("Malformed checkpoint bookmarks")
    if label in bookmarks:
        raise CheckpointError(f"Checkpoint name already exists: {label}")
    base = resolve_checkpoint(registry, parent)
    if model is None:
        if reflections is not None:
            raise CheckpointError("--mtz may be used only when importing --model")
        bookmarks[label] = base["id"]
        save_registry(run, registry)
        return checkpoint_record(base, run)

    source_model = model.expanduser().resolve()
    if not source_model.is_file():
        raise CheckpointError(f"Manual checkpoint model does not exist: {source_model}")
    try:
        _, frozen_report = _run_and_report(run)
        validate_model_phosphate_policy(source_model, frozen_report)
    except PhosphateError as exc:
        raise CheckpointError(f"Manual model phosphate validation failed: {exc}") from exc
    source_reflections: Path | None = None
    if reflections is not None:
        source_reflections = reflections.expanduser().resolve()
        if not source_reflections.is_file():
            raise CheckpointError(
                f"Manual checkpoint reflections do not exist: {source_reflections}"
            )
    checkpoint_id = next_checkpoint_id(registry, "manual")
    destination_dir = run / "AutoRefine" / "Checkpoints" / checkpoint_id
    destination_dir.mkdir(parents=True)
    destination_model = destination_dir / source_model.name
    shutil.copyfile(source_model, destination_model)
    observations = base.get("observations")
    if source_reflections is not None:
        destination_reflections = destination_dir / source_reflections.name
        shutil.copyfile(source_reflections, destination_reflections)
        observations = _file_reference(destination_reflections, run)
    payload = {
        "id": checkpoint_id,
        "parent": base["id"],
        "kind": "manual",
        "status": "REVIEW",
        "usable": True,
        "recipe": "manual-import",
        "label": label,
        "model": _file_reference(destination_model, run),
        "observations": observations,
        "phases": base.get("phases"),
        "restraints": base.get("restraints", []),
        **(
            {
                "terminal_geometry_protection":
                    base["terminal_geometry_protection"]
            }
            if "terminal_geometry_protection" in base
            else {}
        ),
        # The coordination restraints come with the restraint files; their
        # record too, so the same metal cannot be restrained twice (2026-11-24).
        **(
            {"metal_coordination": base["metal_coordination"]}
            if base.get("metal_coordination")
            else {}
        ),
        "metrics": {},
        "compatibility": {
            "validated": False,
            "source": str(source_model),
            "review_required": True,
        },
    }
    append_checkpoint(run, registry, payload, select=False)
    bookmarks = registry.setdefault("bookmarks", {})
    assert isinstance(bookmarks, dict)
    bookmarks[label] = checkpoint_id
    save_registry(run, registry)
    return checkpoint_record(payload, run)


def _inherited_paths(
    item: Mapping[str, object],
    run_directory: Path | None = None,
    *,
    verify_phase_checksum: bool,
) -> dict[str, object]:
    run = run_directory.expanduser().resolve() if run_directory is not None else None
    model = _ref_path(item.get("model"), "model", run)
    observations = _ref_path(item.get("observations"), "observations", run)
    phase = _ref_path(
        item.get("phases"),
        "phase data",
        run,
        required=False,
        verify_checksum=verify_phase_checksum,
    )
    restraint_values = item.get("restraints")
    restraints = [
        _ref_path(value, "restraint", run)
        for value in restraint_values
        if isinstance(restraint_values, Sequence) and isinstance(value, Mapping)
    ] if isinstance(restraint_values, list) else []
    return {
        "model": model,
        "observations": observations,
        "phases": phase,
        "phase_metadata": item.get("phases"),
        "restraints": [path for path in restraints if path is not None],
    }


def inherited_paths(
    item: Mapping[str, object], run_directory: Path | None = None
) -> dict[str, object]:
    """Resolve and verify the frozen files inherited by a refinement child."""
    return _inherited_paths(
        item,
        run_directory,
        verify_phase_checksum=True,
    )


def _inherited_paths_for_autorefine(
    item: Mapping[str, object], run_directory: Path
) -> dict[str, object]:
    """Resolve inputs while deferring phase hashing to AutoRefine's stable check."""
    return _inherited_paths(
        item,
        run_directory,
        verify_phase_checksum=False,
    )


__all__ = [
    "AUTOMATIC_STATUSES",
    "CheckpointError",
    "CheckpointRecord",
    "add_checkpoint",
    "append_checkpoint",
    "checkpoint_record",
    "inherited_paths",
    "initialize_registry",
    "list_checkpoints",
    "next_checkpoint_id",
    "read_registry",
    "registry_path",
    "resolve_checkpoint",
    "save_registry",
    "select_checkpoint",
]
