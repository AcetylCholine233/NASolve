"""Resolve and launch stage-appropriate Coot views in isolated run pens."""

from __future__ import annotations

import json
import math
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
class MetalMerge:
    """Where the in-Coot "merge metal sites and save" command writes."""
    output_directory: Path
    merged_stem: str  # "<starting model>-merged"; ".pdb" (or "-2.pdb", ...) is added
    warnings: tuple[str, ...] = ()


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
    metal_merge: MetalMerge | None = None


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
    metal_merge: MetalMerge | None = None


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


def _site_position_warnings(model: Path, record: Mapping[str, object]) -> tuple[str, ...]:
    """Sites AutoSol placed next to a symmetry copy rather than the model itself.

    The search records each site's closest residue with crystal symmetry
    applied; if the site as placed is clearly farther from the model, merging
    its raw coordinates would put the ion away from the bases it binds.
    """
    atoms = []
    try:
        lines = model.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ()
    for line in lines:
        if line[0:6] in ("ATOM  ", "HETATM") and len(line) >= 54:
            try:
                atoms.append(tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))))
            except ValueError:
                continue
    sites = record.get("sites")
    if not atoms or not isinstance(sites, list):
        return ()
    warnings = []
    for site in sites:
        if not isinstance(site, Mapping):
            continue
        coordinates = site.get("coordinates")
        nearest = site.get("nearest_residues")
        if (
            not isinstance(coordinates, list) or len(coordinates) != 3
            or not isinstance(nearest, list) or not nearest
            or not isinstance(nearest[0], Mapping)
        ):
            continue
        try:
            xyz = tuple(float(value) for value in coordinates)
            with_symmetry = float(nearest[0]["distance"])
        except (KeyError, TypeError, ValueError):
            continue
        as_placed = min(math.dist(atom, xyz) for atom in atoms)
        if as_placed > with_symmetry + 0.5:
            closest = nearest[0]
            warnings.append(
                f"{site.get('element')} site {site.get('index')} is {as_placed:.1f} A from "
                f"the model as placed but {with_symmetry:.1f} A from "
                f"{closest.get('chain')}:{closest.get('residue_number')} in a symmetry "
                "copy; move it next to the model before refining"
            )
    return tuple(warnings)


def _metal_view(
    anchor: Path, record: Mapping[str, object], label: str, output_directory: Path
) -> ViewProfile:
    """Phaser model + AutoSol metal sites over the Phaser map for one search."""
    inputs = record.get("inputs")
    outputs = record.get("outputs")
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise CootViewError("The metal-search record is malformed")
    model = _reported_file(inputs.get("phaser_model"), "Phaser model", anchor)
    map_path = _reported_file(inputs.get("phaser_map"), "Phaser map", anchor)
    heavy_atoms = _reported_file(
        outputs.get("heavy_atom_model"), "AutoSol heavy-atom model", anchor
    )
    return ViewProfile(
        "metal",
        model,
        map_path,
        (heavy_atoms,),
        (),
        f"Phaser solution with AutoSol metal sites ({label})",
        "Phaser MR map",
        metal_merge=MetalMerge(
            output_directory=output_directory,
            merged_stem=f"{model.stem}-merged",
            warnings=_site_position_warnings(model, record),
        ),
    )


def _metal_run_profile(metal_run: Path, report: Mapping[str, object]) -> ViewProfile:
    """A ``DATASET/AutoSolMetal/run_NNN`` directory opened directly."""
    if report.get("status") != "METAL_SITES_FOUND":
        raise CootViewError(
            f"This metal search did not find sites ({report.get('status')}); "
            "there is nothing to show"
        )
    return _metal_view(metal_run, report, f"AutoSolMetal/{metal_run.name}", metal_run)


def _metal_profile(run: Path, report: Mapping[str, object]) -> ViewProfile:
    """The latest metal search with sites that was made from this AutoMR run."""
    block = report.get("metal_autosol")
    searches = block.get("searches") if isinstance(block, Mapping) else None
    found = [
        search for search in (searches if isinstance(searches, list) else [])
        if isinstance(search, Mapping) and search.get("status") == "METAL_SITES_FOUND"
    ]
    if found:
        search = found[-1]
        report_path = _reported_file(search.get("report"), "metal-search report", run)
        expected = search.get("report", {}).get("sha256") if isinstance(
            search.get("report"), Mapping) else None
        _verify_report_checksum(report_path, expected, "metal-search report")
        metal_report = _read_json(report_path, "metal-search report")
        return _metal_run_profile(report_path.parent, metal_report)
    # Schema 1: attempts were recorded (and stored) inside the AutoMR run.
    attempts = block.get("attempts") if isinstance(block, Mapping) else None
    legacy = [
        attempt for attempt in (attempts if isinstance(attempts, list) else [])
        if isinstance(attempt, Mapping) and attempt.get("status") == "METAL_SITES_FOUND"
    ]
    if not legacy:
        raise CootViewError(
            "No metal search from this run found sites; run "
            "'nasolve autosol RUN --metal ELEMENT'"
        )
    attempt = legacy[-1]
    name = str(attempt.get("attempt"))
    attempt_directory = run / "AutoSolMetal" / name
    return _metal_view(
        run, attempt, name, attempt_directory if attempt_directory.is_dir() else run
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
    if report.get("workflow") == "metal_autosol":
        # DATASET/AutoSolMetal/run_NNN opened directly.
        if checkpoint is not None or stage not in {None, "metal"}:
            raise CootViewError("A metal-search directory can only be shown as the metal view")
        return _metal_run_profile(run, report)
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


def _metal_merge_script(model: Path, sites: Path, merge: MetalMerge) -> str:
    """In-Coot command: merge the (possibly edited) sites into a copy of the model.

    Adds a "NASolve" toolbar button with "Merge metal sites and save" and the
    Python function ``nasolve_merge_metal_sites()``. The open molecules are
    used as they are, so sites deleted or moved in Coot stay deleted or moved.
    Occupancies above 1 (AutoSol's relative scale) are set to 1.0. Nothing is
    overwritten: an existing name gets "-2", "-3", ... appended.
    """
    settings = {
        "model": str(model),
        "sites": str(sites),
        "output_directory": str(merge.output_directory),
        "stem": merge.merged_stem,
        "warnings": list(merge.warnings),
    }
    return (
        "\n\n# --- NASolve: merge metal sites into the starting model ---\n"
        f"NASOLVE_MERGE = {json.dumps(settings)}\n"
        '''

def _nasolve_find_molecule(path):
    target = os.path.realpath(path)
    for imol in range(coot.graphics_n_molecules()):
        if coot.is_valid_model_molecule(imol):
            name = coot.molecule_name(imol) or ""
            if os.path.realpath(name) == target:
                return imol
    return -1


def _nasolve_merged_path():
    directory = NASOLVE_MERGE["output_directory"]
    stem = NASOLVE_MERGE["stem"]
    candidate = os.path.join(directory, stem + ".pdb")
    number = 2
    while os.path.exists(candidate):
        candidate = os.path.join(directory, "%s-%d.pdb" % (stem, number))
        number += 1
    return candidate


def _nasolve_cap_occupancies(imol):
    import coot_utils
    changes = []
    notes = []
    n_atoms = 0
    for chain_id in coot_utils.chain_ids(imol):
        for serial in range(coot.chain_n_residues(chain_id, imol)):
            resno = coot.seqnum_from_serial_number(imol, chain_id, serial)
            ins_code = coot.insertion_code_from_serial_number(imol, chain_id, serial)
            for atom in coot.residue_info_py(imol, chain_id, resno, ins_code) or []:
                n_atoms += 1
                name, alt_conf = atom[0][0], atom[0][1]
                occupancy = atom[1][0]
                if occupancy > 1.0:
                    changes.append([imol, chain_id, resno, ins_code, name, alt_conf, "occ", 1.0])
                    notes.append("%s %s%s: occupancy %.2f set to 1.00"
                                 % (name.strip(), chain_id, resno, occupancy))
    if changes:
        coot.set_atom_attributes_py(changes)
    return n_atoms, notes


def nasolve_merge_metal_sites(*_ignored):
    """Merge the open metal sites into a copy of the starting model and save it."""
    imol_model = _nasolve_find_molecule(NASOLVE_MERGE["model"])
    imol_sites = _nasolve_find_molecule(NASOLVE_MERGE["sites"])
    if imol_model < 0 or imol_sites < 0:
        message = ("NASolve: the starting model or overall_best_ha_pdb.pdb is no longer "
                   "open in Coot, so nothing was merged.")
        coot.info_dialog(message)
        print("NASOLVE_MERGED_MODEL " + json.dumps({"status": "missing-molecule"}), flush=True)
        return None
    sites_copy = coot.copy_molecule(imol_sites)
    n_sites, notes = _nasolve_cap_occupancies(sites_copy)
    merged = coot.copy_molecule(imol_model)
    coot.merge_molecules_py([sites_copy], merged)
    coot.close_molecule(sites_copy)
    output = _nasolve_merged_path()
    coot.write_pdb_file(merged, output)
    if not os.path.isfile(output):
        coot.info_dialog("NASolve: Coot did not write " + output)
        print("NASOLVE_MERGED_MODEL " + json.dumps({"status": "write-failed", "path": output}),
              flush=True)
        return None
    coot.set_molecule_name(merged, output)
    lines = ["NASolve saved the merged model:", output,
             "%d metal-site atom(s) merged into a copy of the starting model." % n_sites]
    lines += notes
    if NASOLVE_MERGE["warnings"]:
        lines.append("Check before refining:")
        lines += NASOLVE_MERGE["warnings"]
    coot.info_dialog("\\n".join(lines))
    print("NASOLVE_MERGED_MODEL " + json.dumps({
        "status": "saved", "path": output, "sites": n_sites, "occupancy_changes": notes,
        "warnings": NASOLVE_MERGE["warnings"]}), flush=True)
    return output


def _nasolve_add_merge_button():
    try:
        import coot_gui
        menu = coot_gui.attach_module_menu_button("NASolve")
        coot_gui.add_simple_action_to_menu(
            menu, "Merge metal sites and save", "nasolve_merge_metal_sites",
            lambda _action, _parameter: nasolve_merge_metal_sites())
        print("NASOLVE_MERGE_READY toolbar", flush=True)
    except Exception as exc:
        print("NASOLVE_MERGE_READY scripting-only " + str(exc), flush=True)
    return False


try:
    from gi.repository import GLib
    GLib.timeout_add(500, _nasolve_add_merge_button)
except Exception:
    _nasolve_add_merge_button()
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
    startup_script = _anomalous_map_script(anomalous_paths)
    if profile.metal_merge is not None and profile.extra_model_paths:
        # One startup script: the merge command is appended to the map loader.
        startup_script += _metal_merge_script(
            profile.model_path, profile.extra_model_paths[0], profile.metal_merge
        )
    try:
        map_script.write_text(startup_script, encoding="utf-8")
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
                "metal_merge": (
                    {
                        "output_directory": str(profile.metal_merge.output_directory),
                        "merged_name": profile.metal_merge.merged_stem + ".pdb",
                        "warnings": list(profile.metal_merge.warnings),
                        "log_marker": "NASOLVE_MERGED_MODEL",
                    }
                    if profile.metal_merge is not None else None
                ),
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
        metal_merge=profile.metal_merge,
    )


__all__ = [
    "CootViewError",
    "CootViewResult",
    "MetalMerge",
    "ViewProfile",
    "find_last_run",
    "launch_coot_view",
    "resolve_run",
    "resolve_view_profile",
]
