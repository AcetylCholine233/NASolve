"""Run frozen campaign selections through the existing scientific stage engines.

This module does not schedule work, select review checkpoints, or remember a
workspace. The campaign coordinator owns concurrency, durable receipts, and
interruption policy. A successful Python call is not a scientific acceptance:
the exact stage status is returned for the coordinator to apply its stop gate.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable

from .automr import prepare_automr
from .automr_input import DatasetFiles, ResolvedAutoMRInput, format_intent, normalize_frame
from .autorefine import execute_autorefine
from .autosol import execute_autosol
from .campaign_records import path_in
from .campaigns import CampaignError, _contained, _file_identity, _safe_relative
from .config import load_config
from .coot_runtime import CootDiscoveryError, discover_coot
from .phenix_runtime import PhenixInstallation, discover_phenix
from .phaser import execute_phaser
from .postmr import prepare_postmr
from .residue_aliases import ResolvedLigand


STAGES = ("preflight", "phaser", "postmr", "autosol", "autorefine")


class CampaignStageError(CampaignError):
    """A frozen selection or stage prerequisite is not safe to execute."""


def _object(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise CampaignStageError(f"{label} must be an object")
    return value


def _read_report(path: Path) -> dict[str, Any]:
    try:
        return _object(json.loads(path.read_text(encoding="utf-8")), "Run report")
    except (OSError, ValueError) as exc:
        raise CampaignStageError(f"Could not read run report {path}: {exc}") from exc


def _reference(root: Path, value: object, *, dataset: str | None = None) -> Path:
    reference = _object(value, "Frozen campaign reference")
    if reference.get("anchor") != "campaign":
        raise CampaignStageError("Frozen input has no campaign anchor")
    relative = Path(_safe_relative(reference.get("relative_path"), "Frozen input"))
    allowed = Path("NASolveCampaign/resources") if dataset is None else Path(dataset)
    if not relative.is_relative_to(allowed):
        raise CampaignStageError(f"Frozen input is outside its declared location: {relative}")
    path = _contained(root / allowed, root / relative, "Frozen input")
    if _file_identity(path) != (reference.get("sha256"), reference.get("size")):
        raise CampaignStageError(f"Frozen input checksum or size changed: {relative}")
    return path


def _copy_resource(root: Path, reference: object, destination: Path) -> None:
    source = _reference(root, reference)
    data = source.read_bytes()
    expected = _object(reference, "Frozen resource")
    if (hashlib.sha256(data).hexdigest(), len(data)) != (expected["sha256"], expected["size"]):
        raise CampaignStageError("Frozen resource changed while being copied")
    with destination.open("xb") as handle:
        handle.write(data)


def _ligand(value: object) -> ResolvedLigand:
    fields = _object(value, "Frozen ligand")
    if (not isinstance(fields.get("token"), str)
            or not isinstance(fields.get("ligand_code"), str)
            or not isinstance(fields.get("used_alias"), bool)):
        raise CampaignStageError("Frozen ligand identity is malformed")
    return ResolvedLigand(fields["token"], fields["ligand_code"], fields["used_alias"])


def _pair(value: object) -> tuple[ResolvedLigand, ResolvedLigand]:
    if not isinstance(value, list) or len(value) != 2:
        raise CampaignStageError("Frozen ordered ligand pair is malformed")
    return _ligand(value[0]), _ligand(value[1])


def _frozen_selection(root: Path, dataset: dict[str, Any], attempt: Path) -> ResolvedAutoMRInput:
    """Materialize the chosen model without consulting the original catalogue."""
    if dataset.get("status") != "DISCOVERED":
        raise CampaignStageError("Only a discovered campaign dataset can enter preflight")
    name = dataset["id"]
    effective = _object(dataset.get("effective_config"), "Frozen effective configuration")
    inputs = _object(dataset.get("inputs"), "Frozen inputs")
    if effective.get("mode") != "standard" or effective.get("frame") != "W":
        raise CampaignStageError("Campaign execution supports only frozen standard W selections")
    model_name = Path(_safe_relative(effective.get("model_name"), "Frozen model name"))
    if len(model_name.parts) != 1 or model_name.suffix.lower() != ".pdb":
        raise CampaignStageError("Frozen model name must be one PDB filename")
    if effective.get("model") != inputs.get("model"):
        raise CampaignStageError("Frozen selected model and input reference disagree")
    if effective.get("frame_sequence") != inputs.get("frame_sequence"):
        raise CampaignStageError("Frozen frame sequence references disagree")
    frozen = attempt / "frozen"
    frozen.mkdir()
    model = frozen / model_name
    _copy_resource(root, inputs.get("model"), model)
    frame_sequence_source = None
    if inputs.get("frame_sequence") is not None:
        frame_sequence_source = frozen / "seq_base.txt"
        _copy_resource(root, inputs["frame_sequence"], frame_sequence_source)
    sequence_reference_label = effective.get("sequence_reference")
    sequence_reference_record = inputs.get("sequence_reference")
    if (sequence_reference_label is None) != (sequence_reference_record is None):
        raise CampaignStageError("Frozen sequence-reference selector and resource disagree")
    sequence_reference = None
    if sequence_reference_record is not None:
        if (
            not isinstance(sequence_reference_label, str)
            or not sequence_reference_label.strip()
            or any(char in sequence_reference_label for char in "\r\n\x00")
        ):
            raise CampaignStageError("Frozen sequence-reference selector is malformed")
        sequence_reference = frozen / "sequence_reference.json"
        _copy_resource(root, sequence_reference_record, sequence_reference)
    files = DatasetFiles(
        root=root / name,
        reflections=_reference(root, inputs.get("reflections"), dataset=name),
        metadata=_reference(root, inputs.get("metadata"), dataset=name),
        summary=(
            _reference(root, inputs["summary"], dataset=name)
            if inputs.get("summary") is not None else None
        ),
    )
    resolved = ResolvedAutoMRInput(
        dataset=files, mode="standard", frame=normalize_frame("W"),
        pair_text=effective["pair"], pair=_pair(effective["pair_ligands"]),
        model=model, model_source=effective["model_source"],
        model_pair=(
            _pair(effective["model_pair"])
            if effective.get("model_pair") is not None else None
        ),
        exact_pair_model=effective["exact_pair_model"],
        catalogue_warnings=tuple(effective["catalogue_warnings"]),
        allow_p1_standard=effective["allow_p1_standard"], mirror=effective["mirror"],
        sequences=dict(effective["sequences"]), sequence_file=None,
        sequence_reference=sequence_reference,
        sequence_reference_label=sequence_reference_label,
        sequence_thread=(
            dict(effective["sequence_thread"])
            if effective.get("sequence_thread") is not None else None
        ),
        model_selector=effective.get("model_selector"),
        model_family=effective.get("model_family"),
        model_provider=(
            dict(effective["model_provider"])
            if effective.get("model_provider") is not None else None
        ),
        frame_sequence_source=frame_sequence_source,
        mutations={site: _ligand(ligand) for site, ligand in effective["mutations"].items()},
        config_source=None,
        allow_op3_sites=tuple(effective.get("allow_op3_sites", [])),
        phosphate_intent=effective.get("phosphate_intent"),
        backbone_sites=dict(effective.get("backbones", {})),
        allow_unreviewed_backbone=bool(effective.get("allow_unreviewed_backbone", False)),
    )
    config = frozen / "nasolve.input.txt"
    with config.open("x", encoding="utf-8") as handle:
        handle.write(format_intent(resolved))
    return replace(resolved, config_source=config)


def _program(phenix: PhenixInstallation, name: str) -> Path:
    executable = phenix.executables.get(name)
    if executable is None:
        raise CampaignStageError(f"The discovered Phenix installation has no {name}")
    return executable


def _postmr_report(report: dict[str, Any]) -> dict[str, Any]:
    postmr = _object(report.get("postmr"), "Completed PostMR report")
    if postmr.get("status") != "POSTMR_READY":
        raise CampaignStageError("This stage requires POSTMR_READY")
    return postmr


def _autosol_required(report: dict[str, Any]) -> bool:
    anomalous = _object(_postmr_report(report).get("anomalous"), "PostMR anomalous report")
    candidates = anomalous.get("candidates")
    required = anomalous.get("autosol_required")
    if (not isinstance(candidates, list) or not isinstance(required, bool)
            or required != bool(candidates)
            or not all(isinstance(candidate, dict) for candidate in candidates)):
        raise CampaignStageError("PostMR anomalous eligibility is inconsistent")
    return required


def _relative_file(root: Path, path: Path) -> str:
    resolved = _contained(root, path, "Stage artifact")
    if not resolved.is_file():
        raise CampaignStageError(f"Expected stage artifact is missing: {path}")
    return resolved.relative_to(root).as_posix()


def _tree_files(root: Path, directory: Path) -> list[str]:
    """Receipt every stage-owned regular file; external symlinks are refused."""
    paths: list[str] = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            # A receipt must not silently switch to data outside the stage tree.
            _contained(directory, path, "Stage artifact")
        if path.is_file():
            paths.append(_relative_file(root, path))
    return paths


def _snapshot_report(root: Path, attempt: Path, run: Path) -> str:
    snapshot = attempt / "run-report.json"
    with snapshot.open("xb") as handle:
        handle.write((run / "report.json").read_bytes())
    return _relative_file(root, snapshot)


def execute_stage(
    root: Path,
    dataset: dict[str, Any],
    policy: dict[str, Any],
    stage: str,
    attempt_dir: Path,
    run: Path | None = None,
    *,
    phenix_root: str | None = None,
    on_run_allocated: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    """Execute exactly one stage and return its scientific outcome and artifacts."""
    if stage not in STAGES:
        raise CampaignStageError(f"Unknown campaign stage: {stage}")
    root = root.expanduser().resolve()
    name = Path(_safe_relative(dataset.get("id"), "Campaign dataset"))
    if len(name.parts) != 1 or dataset.get("relative_path") != str(name):
        raise CampaignStageError("Campaign dataset identity is malformed")
    _contained(root, root / name, "Campaign dataset")
    attempt = _contained(root / "NASolveCampaign", attempt_dir, "Campaign attempt")
    if not attempt.is_dir():
        raise CampaignStageError("Campaign attempt directory does not exist")
    if stage == "preflight":
        if run is not None:
            raise CampaignStageError("Preflight cannot replace an existing campaign run")
        # Check before prepare_automr allocates its numbered run. Its standalone
        # allocator intentionally accepts ordinary paths; a campaign must not
        # follow a substituted AutoMR directory into another dataset or tree.
        path_in(root, name.as_posix() + "/AutoMR")
    else:
        if run is None:
            raise CampaignStageError(f"{stage} requires an allocated campaign run")
        run = _contained(root / name / "AutoMR", run, "Campaign run")
        if run.parent != root / name / "AutoMR" or not re.fullmatch(r"run_\d{3,}", run.name):
            raise CampaignStageError("Campaign run must be one numbered AutoMR directory")
    config = load_config()
    report = _read_report(run / "report.json") if run is not None else None
    if report is not None and report.get("workflow") != "automr":
        raise CampaignStageError("Campaign stage requires an AutoMR run")
    if stage == "autosol" and not _autosol_required(report):
        snapshot = _snapshot_report(root, attempt, run)
        return {
            "status": "SKIPPED", "message": "PostMR found no nucleotide heavy atom; AutoSol is not required",
            "run": run.relative_to(root).as_posix(), "artifacts": [snapshot],
            "run_report_snapshot": snapshot, "tool_versions": {},
        }
    if stage == "autorefine" and _autosol_required(report):
        autosol = _object(report.get("autosol"), "Accepted AutoSol report")
        if (autosol.get("status") != "AUTOSOL_READY"
                or autosol.get("use_for_refinement") is not True):
            raise CampaignStageError("The anomalous branch requires AUTOSOL_READY before refinement")
    phenix = discover_phenix(config, explicit=phenix_root)
    versions = {"phenix": phenix.version}
    extra: dict[str, Any] = {}
    if stage == "preflight":
        resolved = _frozen_selection(root, dataset, attempt)
        result = prepare_automr(
            resolved.dataset.root, resolved_input=resolved,
            mtz_dump_executable=_program(phenix, "phenix.mtz.dump"),
            phenix_environment=phenix.environment,
            on_run_allocated=on_run_allocated,
        )
        run = result.run_directory
        stage_directory = run / "Model"
        required = [run / "nasolve.input.txt", stage_directory / "input_model.pdb",
                    stage_directory / "assessment.json"]
        extra["artifacts"] = _tree_files(root, attempt / "frozen")
    elif stage == "phaser":
        result = execute_phaser(
            run / "report.json", _program(phenix, "phenix.phaser"),
            environment=phenix.environment, phenix_version=phenix.version,
        )
        stage_directory = result.phaser_directory
        required = [result.parameter_path, result.log_path]
        if result.status in {"MR_SUCCESS", "MR_REVIEW"}:
            required.extend([stage_directory / "mr_solution.pdb", stage_directory / "mr_solution.mtz"])
        extra.update(tfz=result.tfz, llg=result.llg)
    elif stage == "postmr":
        if report.get("status") != "MR_SUCCESS":
            raise CampaignStageError("Campaign PostMR requires MR_SUCCESS; review is an inspection stop")
        coot = None
        try:
            coot = discover_coot(config)
        except CootDiscoveryError:
            pass  # The engine requires Coot only if its mutation plan needs it.
        if coot is not None:
            versions["coot"] = coot.version
        postmr_policy = _object(policy.get("postmr"), "PostMR policy")
        if set(postmr_policy) != {"modified_pairs_only"} or type(postmr_policy["modified_pairs_only"]) is not bool:
            raise CampaignStageError("Unsupported campaign PostMR policy")
        result = prepare_postmr(
            run, _program(phenix, "phenix.ready_set"),
            coot_executable=coot.executable if coot else None,
            environment=phenix.environment,
            modified_pairs_only=postmr_policy["modified_pairs_only"],
        )
        stage_directory = result.postmr_directory
        required = [result.model_path, result.readyset_log, *result.restraint_paths]
    elif stage == "autosol":
        autosol_policy = _object(policy.get("autosol"), "AutoSol policy")
        if autosol_policy != {"policy": "when-anomalous", "on_unaccepted": "inspect"}:
            raise CampaignStageError("Unsupported campaign AutoSol policy")
        result = execute_autosol(
            run, _program(phenix, "phenix.autosol"), _program(phenix, "phenix.mtz.dump"),
            environment=phenix.environment,
        )
        stage_directory = result.autosol_directory
        required = [result.log_path] if result.log_path.is_file() else []
        if result.status == "AUTOSOL_WARNING":
            warning = _read_report(result.report_path)
            extra["message"] = (
                f"AutoSol phases were not accepted: {warning.get('failure_reason', result.message)}; "
                "campaign stopped for inspection"
            )
        elif not required:
            raise CampaignStageError("Completed AutoSol stage has no console log")
        if result.status == "AUTOSOL_READY":
            if result.heavy_atom_model is None or result.refinement_data is None:
                raise CampaignStageError("AUTOSOL_READY has no heavy-atom model or refinement data")
            required.extend([result.heavy_atom_model, result.refinement_data])
        extra["matched_distance"] = result.matched_distance
    else:
        refine_policy = _object(policy.get("autorefine"), "AutoRefine policy")
        if refine_policy != {"recipe": "AutoRefine/default", "cycles": 5}:
            raise CampaignStageError("Unsupported campaign AutoRefine recipe")
        if (run / "AutoRefine").exists():
            raise CampaignStageError("Campaign first refinement refuses an existing AutoRefine directory")
        result = execute_autorefine(
            run, _program(phenix, "phenix.refine"), _program(phenix, "phenix.mtz.dump"),
            phenix_version=phenix.version, environment=phenix.environment,
            from_checkpoint="postmr", recipe=refine_policy["recipe"],
            macro_cycles=refine_policy["cycles"],
        )
        stage_directory = result.round_directory
        required = [result.log_path]
        if result.status in {"AUTOREFINE_READY", "AUTOREFINE_REVIEW", "AUTOREFINE_ANOMALOUS_FALLBACK"}:
            if result.model_path is None or result.map_coefficients is None:
                raise CampaignStageError("Completed refinement has no model or map coefficients")
            required.extend([result.model_path, result.map_coefficients])
        registry_snapshot = attempt / "checkpoints.json"
        with registry_snapshot.open("xb") as handle:
            handle.write((run / "AutoRefine" / "checkpoints.json").read_bytes())
        required.append(registry_snapshot)
        extra.update(checkpoint=result.checkpoint_id, selected_as_current=result.selected_as_current,
                     statistics=result.statistics)
    required.append(result.report_path)
    # report.json changes in later stages. Preserve its exact boundary contents
    # for the coordinator rather than treating the live report as immutable.
    snapshot = _snapshot_report(root, attempt, run)
    artifacts = _tree_files(root, stage_directory)
    artifacts.extend(extra.pop("artifacts", []))
    artifacts.extend(_relative_file(root, path) for path in required if path != run / "report.json")
    artifacts.append(snapshot)
    return {
        "status": result.status, "message": result.message,
        "run": run.relative_to(root).as_posix(),
        "artifacts": sorted(set(artifacts)), "run_report_snapshot": snapshot,
        "tool_versions": versions, **extra,
    }
