"""Resolve and launch stage-appropriate Coot views in isolated run pens."""

from __future__ import annotations

import json
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping

from .checkpoints import READABLE_SCHEMA_VERSIONS
from .model_assessment import file_sha256
from .ligand_profiles import is_modification_only
from .run_context import resolve_artifact_path


class CootViewError(RuntimeError):
    """Raised when a requested result cannot be opened safely."""


_CHECKPOINT_COMPONENT = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")


@dataclass(frozen=True)
class ViewProfile:
    stage: str
    model_path: Path
    map_path: Path
    extra_model_paths: tuple[Path, ...]
    dictionary_paths: tuple[Path, ...]
    source: str
    map_source: str
    checkpoint_id: str | None = None
    anomalous_map_paths: tuple[Path, ...] = ()


@dataclass(frozen=True)
class CootViewResult:
    run_directory: Path
    stage: str
    working_directory: Path
    model_path: Path
    map_path: Path
    extra_model_paths: tuple[Path, ...]
    dictionary_paths: tuple[Path, ...]
    log_path: Path
    command: tuple[str, ...]
    pid: int
    source: str
    map_source: str
    checkpoint_id: str | None


def _required_file(path: Path, description: str) -> Path:
    if not path.is_file():
        raise CootViewError(f"{description} is missing: {path}")
    return path.resolve()


def _checkpoint_component(value: str) -> str:
    """Return a checkpoint ID that is safe to use as one directory name."""
    if _CHECKPOINT_COMPONENT.fullmatch(value) is None:
        raise CootViewError(
            "Checkpoint IDs used for Coot must be a single safe path component"
        )
    return value


def _reject_symlink_components(root: Path, target: Path) -> None:
    """Reject an existing symlink anywhere beneath a machine-local output path."""
    try:
        relative = target.relative_to(root)
    except ValueError as exc:
        raise CootViewError(f"Coot working path escapes the selected run: {target}") from exc
    candidate = root
    for part in relative.parts:
        candidate = candidate / part
        if candidate.is_symlink():
            raise CootViewError(
                f"Coot working path contains a symbolic link: {candidate}"
            )


def _require_scratch_containment(scratch_root: Path, path: Path) -> None:
    """Verify that a created scratch path still resolves beneath its lexical root."""
    try:
        path.resolve().relative_to(scratch_root)
    except ValueError as exc:
        raise CootViewError(f"Coot working path escapes its scratch directory: {path}") from exc


def _reported_file(value: object, description: str, run: Path) -> Path:
    if isinstance(value, Mapping):
        _validate_report_checksum(value.get("sha256"), description)
    path = resolve_artifact_path(value, run)
    if path is None:
        raise CootViewError(
            f"Reported {description} is missing or failed checksum validation"
        )
    return _required_file(path, description)


def _validate_report_checksum(expected: object, description: str) -> None:
    if expected is None:
        return
    if (
        not isinstance(expected, str)
        or re.fullmatch(r"[0-9a-fA-F]{64}", expected) is None
    ):
        raise CootViewError(f"Reported {description} checksum is malformed")


def _verify_report_checksum(path: Path, expected: object, description: str) -> None:
    _validate_report_checksum(expected, description)
    if isinstance(expected, str) and file_sha256(path) != expected.casefold():
        raise CootViewError(f"Reported {description} failed checksum validation")


def _require_checkpoint_checksum(value: object, description: str) -> None:
    if not isinstance(value, Mapping) or not isinstance(value.get("anchor"), str):
        raise CootViewError(f"Schema-2 checkpoint {description} is not portable")
    expected = value.get("sha256")
    if (
        not isinstance(expected, str)
        or re.fullmatch(r"[0-9a-fA-F]{64}", expected) is None
    ):
        raise CootViewError(f"Schema-2 checkpoint {description} checksum is malformed")


def _read_json(path: Path, description: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CootViewError(f"Could not read {description} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise CootViewError(f"Malformed {description}: {path}")
    return value


def find_last_run(dataset: Path) -> Path:
    """Return the highest numbered run with a readable NASolve report."""
    root = dataset.expanduser().resolve()
    automr = root / "AutoMR"
    candidates: list[tuple[int, Path]] = []
    if automr.is_dir():
        for path in automr.iterdir():
            match = re.fullmatch(r"run_(\d+)", path.name)
            report_path = path / "report.json"
            if not match or not path.is_dir() or not report_path.is_file():
                continue
            try:
                report = json.loads(report_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if isinstance(report, Mapping) and report.get("workflow") == "automr":
                candidates.append((int(match.group(1)), path.resolve()))
    if not candidates:
        raise CootViewError(f"No completed NASolve run was found beneath {root}")
    return max(candidates, key=lambda item: item[0])[1]


def resolve_run(target: Path | str, dataset: Path | None = None) -> Path:
    if str(target).casefold() == "last":
        if dataset is None:
            raise CootViewError("'nasolve show last' requires a dataset directory")
        return find_last_run(dataset)
    if dataset is not None:
        raise CootViewError("A dataset is accepted only with 'nasolve show last DATASET'")
    run = Path(target).expanduser().resolve()
    _required_file(run / "report.json", "AutoMR report")
    return run


def _postmr_dictionaries(
    run: Path, report: Mapping[str, object]
) -> tuple[Path, ...]:
    postmr = report.get("postmr")
    if isinstance(postmr, Mapping) and "view_dictionaries" in postmr:
        references = postmr["view_dictionaries"]
        if not isinstance(references, list) or not references:
            raise CootViewError("Malformed frozen PostMR dictionary list")
        dictionaries = []
        for ref in references:
            _require_checkpoint_checksum(ref, "PostMR ligand dictionary")
            dictionaries.append(_reported_file(ref, "PostMR ligand dictionary", run))
        return tuple(dictionaries)
    readyset = postmr.get("readyset") if isinstance(postmr, Mapping) else None
    generated = readyset.get("generated_ligand_cif") if isinstance(readyset, Mapping) else None
    if isinstance(generated, Mapping):
        _validate_report_checksum(generated.get("sha256"), "ReadySet ligand dictionary")
    generated_path = resolve_artifact_path(generated, run)
    if generated_path is not None:
        return (generated_path,)
    if isinstance(generated, Mapping) and generated.get("sha256") is not None:
        raise CootViewError(
            "ReadySet ligand dictionary is missing or failed checksum validation"
        )
    restraints = postmr.get("restraints") if isinstance(postmr, Mapping) else None
    if not isinstance(restraints, list):
        restraints = []
    dictionaries: list[Path] = []
    for value in restraints:
        filename = (
            value.get("relative_path", value.get("absolute_path", value.get("path")))
            if isinstance(value, Mapping) else value
        )
        if not isinstance(filename, str) or not filename.casefold().endswith(".cif"):
            continue
        path = resolve_artifact_path(value, run)
        if path is None:
            raise CootViewError(
                "PostMR ligand dictionary is missing or could not be rebased"
            )
        dictionaries.append(path)
    if generated is not None and not dictionaries:
        raise CootViewError(
            "ReadySet ligand dictionary is missing and no input CIF fallback is available"
        )
    return tuple(dictionaries)


def _phaser_map(run: Path, report: Mapping[str, object]) -> Path:
    execution = report.get("execution")
    phaser = execution.get("phaser") if isinstance(execution, Mapping) else None
    value = phaser.get("solution_mtz") if isinstance(phaser, Mapping) else None
    return (
        _reported_file(value, "Phaser reflections", run)
        if value is not None
        else _required_file(run / "Phaser" / "mr_solution.mtz", "Phaser reflections")
    )


def _phaser_profile(run: Path, report: Mapping[str, object]) -> ViewProfile:
    execution = report.get("execution")
    phaser = execution.get("phaser") if isinstance(execution, Mapping) else None
    model_value = phaser.get("solution_pdb") if isinstance(phaser, Mapping) else None
    model = (
        _reported_file(model_value, "Phaser model", run)
        if model_value is not None
        else _required_file(run / "Phaser" / "mr_solution.pdb", "Phaser model")
    )
    map_path = _phaser_map(run, report)
    return ViewProfile(
        "automr", model, map_path, (), (), "Phaser solution", "Phaser MR map"
    )


def _metal_profile(run: Path, report: Mapping[str, object]) -> ViewProfile:
    """Phaser model + AutoSol metal sites over the Phaser map (latest found attempt)."""
    block = report.get("metal_autosol")
    attempts = block.get("attempts") if isinstance(block, Mapping) else None
    found = [
        attempt for attempt in (attempts if isinstance(attempts, list) else [])
        if isinstance(attempt, Mapping) and attempt.get("status") == "METAL_SITES_FOUND"
    ]
    if not found:
        raise CootViewError(
            "No metal search in this run found sites; run 'nasolve autosol RUN --metal ELEMENT'"
        )
    attempt = found[-1]
    inputs = attempt.get("inputs")
    outputs = attempt.get("outputs")
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise CootViewError("The metal-search record is malformed")
    model = _reported_file(inputs.get("phaser_model"), "Phaser model", run)
    map_path = _reported_file(inputs.get("phaser_map"), "Phaser map", run)
    heavy_atoms = _reported_file(
        outputs.get("heavy_atom_model"), "AutoSol heavy-atom model", run
    )
    return ViewProfile(
        "metal",
        model,
        map_path,
        (heavy_atoms,),
        (),
        f"Phaser solution with AutoSol metal sites ({attempt.get('attempt')})",
        "Phaser MR map",
    )


def _prepared_model(run: Path, report: Mapping[str, object]) -> Path:
    postmr = report.get("postmr")
    if not isinstance(postmr, Mapping) or postmr.get("status") != "POSTMR_READY":
        raise CootViewError("This run has no completed PostMR result")
    model = _reported_file(postmr.get("prepared_model"), "ReadySet model", run)
    _verify_report_checksum(model, postmr.get("prepared_sha256"), "ReadySet model")
    return model


def _accepted_autosol(report: Mapping[str, object]) -> bool:
    autosol = report.get("autosol")
    return (
        isinstance(autosol, Mapping)
        and autosol.get("status") == "AUTOSOL_READY"
        and autosol.get("use_for_refinement") is True
    )


def _postmr_map(run: Path, report: Mapping[str, object]) -> tuple[Path, str]:
    if _accepted_autosol(report):
        density_map = _autosol_map(run, report, allow_missing=True)
        if density_map is not None:
            return density_map, "accepted AutoSol density-modified map"
        reason = "accepted AutoSol has no density-modified map"
    elif isinstance(report.get("autosol"), Mapping):
        reason = "AutoSol phases are not accepted"
    else:
        reason = "no AutoSol result"
    return _phaser_map(run, report), f"Phaser MR map ({reason})"


def _postmr_profile(run: Path, report: Mapping[str, object]) -> ViewProfile:
    model = _prepared_model(run, report)
    map_path, map_source = _postmr_map(run, report)
    return ViewProfile(
        "postmr",
        model,
        map_path,
        (),
        _postmr_dictionaries(run, report),
        "PostMR ReadySet prepared model",
        map_source,
    )


def _autosol_map(
    run: Path, report: Mapping[str, object], *, allow_missing: bool = False
) -> Path | None:
    autosol = report.get("autosol")
    outputs = autosol.get("outputs") if isinstance(autosol, Mapping) else None
    value = outputs.get("density_modified_map") if isinstance(outputs, Mapping) else None
    if value is not None:
        # A declared map is authoritative. Never replace a missing/corrupt one
        # with an unreported product, or use the HL-only refinement data as a map.
        path = _reported_file(value, "AutoSol density-modified map", run)
        _verify_report_checksum(
            path, outputs.get("density_modified_map_sha256"),
            "AutoSol density-modified map",
        )
        return path
    status = outputs.get("density_modified_map_status") if isinstance(outputs, Mapping) else None
    if isinstance(status, str) and status in {"missing", "unavailable"}:
        if allow_missing:
            return None
        raise CootViewError(f"AutoSol density-modified map is {status}")
    if status is not None:
        raise CootViewError(
            f"AutoSol density-modified map has no usable declared source ({status})"
        )
    preferred = [
        path for path in (run / "AutoSol").rglob("overall_best_denmod_map_coeffs.mtz")
        if path.is_file() and "PDS" not in path.parts and "TEMP" not in path.parts
    ]
    if not preferred:
        preferred = [
            path for path in (run / "AutoSol").rglob("resolve_histograms*.mtz")
            if path.is_file() and "PDS" not in path.parts and "TEMP" not in path.parts
        ]
    if not preferred and allow_missing:
        return None
    if len(preferred) != 1:
        raise CootViewError(
            f"Expected one AutoSol density-modified map MTZ, found {len(preferred)}"
        )
    return preferred[0].resolve()


def _autosol_profile(run: Path, report: Mapping[str, object]) -> ViewProfile:
    autosol = report.get("autosol")
    outputs = autosol.get("outputs") if isinstance(autosol, Mapping) else None
    if not isinstance(outputs, Mapping):
        raise CootViewError("This run has no completed AutoSol output")
    heavy_atom = _reported_file(
        outputs.get("heavy_atom_model"), "AutoSol HA model", run
    )
    _verify_report_checksum(
        heavy_atom,
        outputs.get("heavy_atom_model_sha256"),
        "AutoSol HA model",
    )
    model = _prepared_model(run, report)
    map_path = _autosol_map(run, report)
    assert map_path is not None
    status = str(autosol.get("status", "legacy status unrecorded"))
    acceptance = "accepted" if _accepted_autosol(report) else "unaccepted; inspection only"
    return ViewProfile(
        "autosol",
        model,
        map_path,
        (heavy_atom,),
        _postmr_dictionaries(run, report),
        f"PostMR ReadySet prepared model with separate AutoSol HA sites ({status})",
        f"AutoSol density-modified map ({acceptance}; {status})",
    )


def _refinement_map(
    run: Path, checkpoint: Mapping[str, object], schema_version: int
) -> Path:
    outputs = checkpoint.get("outputs")
    map_value = outputs.get("map_coefficients") if isinstance(outputs, Mapping) else None
    if schema_version == 2 and map_value is not None:
        _require_checkpoint_checksum(map_value, "map coefficients")
    reported_map = resolve_artifact_path(map_value, run)
    if reported_map is not None:
        return reported_map
    if map_value is not None:
        raise CootViewError(
            "Checkpoint map coefficients are missing or failed checksum validation"
        )
    if schema_version == 2:
        reflection_value = (
            outputs.get("refinement_reflections")
            if isinstance(outputs, Mapping)
            else None
        )
        _require_checkpoint_checksum(reflection_value, "refinement reflections")
        return _reported_file(reflection_value, "checkpoint refinement reflections", run)
    # Phenix 1.20.x stores 2mFo-DFc/mFo-DFc coefficients in the primary
    # refinement MTZ rather than a separately named *_map_coeffs.mtz file.
    model = _reported_file(checkpoint.get("model"), "checkpoint model", run)
    return _required_file(model.with_suffix(".mtz"), "refinement map coefficients")


def _checkpoint_dictionaries(
    run: Path, checkpoint: Mapping[str, object], schema_version: int
) -> tuple[Path, ...]:
    restraints = checkpoint.get("restraints")
    dictionaries: list[Path] = []
    if isinstance(restraints, list):
        for value in restraints:
            if schema_version == 2:
                _require_checkpoint_checksum(value, "restraint")
            path = resolve_artifact_path(value, run)
            if path is not None:
                if path.suffix.casefold() == ".cif" and not is_modification_only(path):
                    dictionaries.append(path)
                continue
            reported_value = (
                value.get("relative_path", value.get("absolute_path", value.get("path")))
                if isinstance(value, Mapping) else value
            )
            if isinstance(reported_value, str) and reported_value.casefold().endswith(".cif"):
                raise CootViewError(
                    "Checkpoint ligand dictionary is missing or failed checksum validation"
                )
    return tuple(dictionaries)


def _refinement_anomalous_maps(
    run: Path, checkpoint: Mapping[str, object], schema_version: int, map_path: Path
) -> tuple[Path, ...]:
    """Check only this checkpoint's maps and explicitly recorded refinement MTZ."""
    outputs = checkpoint.get("outputs")
    value = outputs.get("refinement_reflections") if isinstance(outputs, Mapping) else None
    if value is None:
        return (map_path,)
    if schema_version == 2:
        _require_checkpoint_checksum(value, "refinement reflections")
    reflections = _reported_file(value, "checkpoint refinement reflections", run)
    return tuple(dict.fromkeys((map_path, reflections)))


def _same_observations(
    run: Path, child: Mapping[str, object], parent: Mapping[str, object], schema_version: int
) -> None:
    paths = []
    for item in (child, parent):
        value = item.get("observations")
        if schema_version == 2:
            _require_checkpoint_checksum(value, "observations")
        paths.append(_reported_file(value, "checkpoint observations", run))
    if paths[0] != paths[1] and file_sha256(paths[0]) != file_sha256(paths[1]):
        raise CootViewError(
            "Manual checkpoint observations differ from its ancestor; no compatible "
            "map coefficients are recorded. Refine this checkpoint before viewing it."
        )


def _checkpoint_map(
    run: Path, report: Mapping[str, object], checkpoint: Mapping[str, object],
    checkpoints: Mapping[str, Mapping[str, object]], schema_version: int,
    visited: frozenset[str] = frozenset(),
) -> tuple[Path, str, tuple[Path, ...]]:
    identifier = str(checkpoint["id"])
    if identifier in visited:
        raise CootViewError("Checkpoint map lineage contains a cycle")
    kind = checkpoint.get("kind")
    if kind == "refinement":
        map_path = _refinement_map(run, checkpoint, schema_version)
        return (
            map_path,
            f"refinement checkpoint {identifier}",
            _refinement_anomalous_maps(run, checkpoint, schema_version, map_path),
        )
    if kind == "postmr":
        map_path, map_source = _postmr_map(run, report)
        return map_path, map_source, (map_path,)
    if kind != "manual":
        raise CootViewError(f"Checkpoint {identifier} has unsupported kind: {kind}")
    parent_id = checkpoint.get("parent")
    parent = checkpoints.get(parent_id) if isinstance(parent_id, str) else None
    if parent is None:
        raise CootViewError(f"Manual checkpoint {identifier} has no recorded map ancestor")
    _same_observations(run, checkpoint, parent, schema_version)
    map_path, map_source, anomalous_paths = _checkpoint_map(
        run, report, parent, checkpoints, schema_version, visited | {identifier}
    )
    return (
        map_path,
        f"inherited {map_source} (not recalculated for manual model)",
        anomalous_paths,
    )


def _autorefine_profile(
    run: Path,
    report: Mapping[str, object],
    checkpoint_id: str | None = None,
) -> ViewProfile:
    registry = _read_json(run / "AutoRefine" / "checkpoints.json", "checkpoint registry")
    schema_version = registry.get("schema_version")
    if type(schema_version) is not int or schema_version not in READABLE_SCHEMA_VERSIONS:
        raise CootViewError(f"Unsupported checkpoint schema: {schema_version}")
    current = checkpoint_id or registry.get("current")
    items = registry.get("checkpoints")
    if not isinstance(current, str) or not isinstance(items, list):
        raise CootViewError("Checkpoint registry has no valid current pointer")
    _checkpoint_component(current)
    checkpoints: dict[str, Mapping[str, object]] = {}
    for item in items:
        if not isinstance(item, Mapping) or not isinstance(item.get("id"), str):
            raise CootViewError("Checkpoint registry contains a malformed checkpoint")
        identifier = item["id"]
        if identifier in checkpoints:
            raise CootViewError(f"Checkpoint registry contains duplicate ID: {identifier}")
        checkpoints[identifier] = item
    checkpoint = checkpoints.get(current)
    if checkpoint is None:
        requested = "requested checkpoint" if checkpoint_id else "current pointer"
        raise CootViewError(f"Checkpoint registry {requested} is unknown: {current}")
    model_ref = checkpoint.get("model")
    if schema_version == 2:
        _require_checkpoint_checksum(model_ref, "model")
    model = _reported_file(model_ref, "checkpoint model", run)
    map_path, map_source, anomalous_paths = _checkpoint_map(
        run, report, checkpoint, checkpoints, schema_version
    )
    selection = "requested" if checkpoint_id is not None else "current"
    status = checkpoint.get("status", "legacy status unrecorded")
    return ViewProfile(
        "autorefine",
        model,
        map_path,
        (),
        _checkpoint_dictionaries(run, checkpoint, schema_version),
        f"{selection} {checkpoint.get('kind')} checkpoint {current} ({status})",
        map_source,
        current,
        anomalous_paths,
    )


def resolve_view_profile(
    run_directory: Path,
    *,
    stage: str | None = None,
    checkpoint: str | None = None,
) -> ViewProfile:
    run = run_directory.expanduser().resolve()
    if checkpoint is not None:
        checkpoint = _checkpoint_component(checkpoint)
    report = _read_json(_required_file(run / "report.json", "AutoMR report"), "run report")
    resolvers = {
        "automr": _phaser_profile,
        "metal": _metal_profile,
        "postmr": _postmr_profile,
        "autosol": _autosol_profile,
        "autorefine": _autorefine_profile,
    }
    if checkpoint is not None:
        if stage not in {None, "autorefine"}:
            raise CootViewError("--checkpoint may be used only with the AutoRefine view")
        return _autorefine_profile(run, report, checkpoint)
    if stage is not None:
        try:
            return resolvers[stage](run, report)
        except KeyError as exc:
            raise CootViewError(f"Unknown viewing stage: {stage}") from exc

    registry_path = run / "AutoRefine" / "checkpoints.json"
    if registry_path.is_file():
        return _autorefine_profile(run, report)
    if _accepted_autosol(report) and _autosol_map(run, report, allow_missing=True) is not None:
        return _autosol_profile(run, report)
    postmr = report.get("postmr")
    if isinstance(postmr, Mapping) and postmr.get("status") == "POSTMR_READY":
        return _postmr_profile(run, report)
    execution = report.get("execution")
    phaser = execution.get("phaser") if isinstance(execution, Mapping) else None
    if isinstance(phaser, Mapping) or (run / "Phaser" / "mr_solution.pdb").is_file():
        return _phaser_profile(run, report)
    raise CootViewError("No completed stage in this run can be opened in Coot")


def _anomalous_map_script(map_paths: tuple[Path, ...]) -> str:
    """Use Coot's MTZ reader in the GUI process; no host crystallography dependency."""
    return (
        "import coot\nimport json\nimport os\n\n"
        f"map_paths = {json.dumps([str(path) for path in map_paths])}\n"
        '''
def report_anomalous_map(status, **details):
    details["status"] = status
    details["labels"] = ["ANOM", "PHANOM"]
    print("NASOLVE_ANOMALOUS_MAP " + json.dumps(details, sort_keys=True), flush=True)


def existing_anomalous_map(path):
    # --auto or the user's Coot configuration may already have opened it.
    for imol in range(coot.graphics_n_molecules()):
        source = coot.mtz_hklin_for_map(imol)
        if source and os.path.realpath(source) == os.path.realpath(path):
            amplitude = coot.mtz_fp_for_map(imol) or ""
            phase = coot.mtz_phi_for_map(imol) or ""
            if (amplitude.rsplit("/", 1)[-1] == "ANOM"
                    and phase.rsplit("/", 1)[-1] == "PHANOM"
                    and not coot.mtz_use_weight_for_map(imol)):
                return imol
    return -1


def load_anomalous_map():
    for path in map_paths:
        if coot.valid_labels(path, "ANOM", "PHANOM", "", 0) != 1:
            continue
        refinement_map = coot.imol_refinement_map()
        scroll_map = coot.scroll_wheel_map()
        try:
            imol = existing_anomalous_map(path)
            reused = imol >= 0
            if not reused:
                # PHANOM already contains the map phase: no extra weights or shift.
                imol = coot.make_and_draw_map(path, "ANOM", "PHANOM", "", 0, 1)
            if imol < 0:
                raise RuntimeError("Coot could not construct ANOM/PHANOM from " + path)
            coot.set_map_is_difference_map(imol, 1)
            coot.set_contour_level_in_sigma(imol, 3.0)
            coot.set_molecule_name(imol, "Anomalous difference (ANOM/PHANOM)")
        finally:
            # An inspection overlay must not change real-space refinement/scrolling.
            if refinement_map >= 0:
                coot.set_imol_refinement_map(refinement_map)
            if scroll_map >= 0:
                coot.set_scroll_wheel_map(scroll_map)
        report_anomalous_map("reused" if reused else "loaded", path=path,
                             molecule=imol, contour_sigma=3.0)
        return
    report_anomalous_map("absent", paths=map_paths)


try:
    load_anomalous_map()
except Exception as exc:
    # Ordinary --auto maps remain available; make optional-map failures visible.
    report_anomalous_map("error", message=str(exc), paths=map_paths)
'''
    )


def launch_coot_view(
    run_directory: Path,
    coot_executable: Path,
    *,
    stage: str | None = None,
    checkpoint: str | None = None,
    environment: Mapping[str, str] | None = None,
    launcher: Callable[..., object] = subprocess.Popen,
) -> CootViewResult:
    """Open a resolved result without allowing Coot to litter the repository."""
    run = run_directory.expanduser().resolve()
    if checkpoint is not None:
        checkpoint = _checkpoint_component(checkpoint)
    profile = resolve_view_profile(run, stage=stage, checkpoint=checkpoint)
    scratch_root = run / "CootGUI"
    working = scratch_root / profile.stage
    if profile.checkpoint_id is not None:
        working = working / profile.checkpoint_id
    backups = working / "backups"
    _reject_symlink_components(run, backups)
    try:
        backups.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise CootViewError(f"Could not create Coot working directory: {exc}") from exc
    _reject_symlink_components(run, backups)
    _require_scratch_containment(scratch_root, working)
    _require_scratch_containment(scratch_root, backups)
    log = working / "coot_gui.log"
    launch_record = working / "launch.json"
    map_script = working / "load_anomalous_map.py"
    for output in (log, launch_record, map_script):
        if output.is_symlink():
            raise CootViewError(f"Refusing to write Coot output through a symlink: {output}")
    command = [
        str(coot_executable.expanduser().resolve()),
        "--pdb",
        str(profile.model_path),
    ]
    for extra in profile.extra_model_paths:
        command.extend(["--pdb", str(extra)])
    for dictionary in profile.dictionary_paths:
        command.extend(["--dictionary", str(dictionary)])
    command.extend(["--auto", str(profile.map_path)])
    anomalous_paths = profile.anomalous_map_paths or (profile.map_path,)
    try:
        map_script.write_text(_anomalous_map_script(anomalous_paths), encoding="utf-8")
    except OSError as exc:
        raise CootViewError(f"Could not write Coot map startup script: {exc}") from exc
    command.extend(["--script", str(map_script)])

    env = dict(os.environ if environment is None else environment)
    env["COOT_BACKUP_DIR"] = str(backups)
    try:
        with log.open("a", encoding="utf-8") as stream:
            process = launcher(
                command,
                cwd=working,
                env=env,
                stdout=stream,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
    except OSError as exc:
        raise CootViewError(f"Could not launch graphical Coot: {exc}") from exc
    pid = getattr(process, "pid", None)
    if not isinstance(pid, int):
        raise CootViewError("Graphical Coot launcher returned no process ID")
    launch_record.write_text(
        json.dumps(
            {
                "pid": pid,
                "stage": profile.stage,
                "run_directory": str(run),
                "checkpoint_id": profile.checkpoint_id,
                "source": profile.source,
                "map_source": profile.map_source,
                "model_path": str(profile.model_path),
                "map_path": str(profile.map_path),
                "anomalous_map": {
                    "status": "check-in-coot",
                    "paths": [str(path) for path in anomalous_paths],
                    "labels": ["ANOM", "PHANOM"],
                    "contour_sigma": 3.0,
                    "script": str(map_script),
                    "log_marker": "NASOLVE_ANOMALOUS_MAP",
                },
                "command": command,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    return CootViewResult(
        run_directory=run,
        stage=profile.stage,
        working_directory=working,
        model_path=profile.model_path,
        map_path=profile.map_path,
        extra_model_paths=profile.extra_model_paths,
        dictionary_paths=profile.dictionary_paths,
        log_path=log,
        command=tuple(command),
        pid=pid,
        source=profile.source,
        map_source=profile.map_source,
        checkpoint_id=profile.checkpoint_id,
    )


__all__ = [
    "CootViewError",
    "CootViewResult",
    "ViewProfile",
    "find_last_run",
    "launch_coot_view",
    "resolve_run",
    "resolve_view_profile",
]
