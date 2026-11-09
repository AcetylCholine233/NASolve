"""AutoMR preflight orchestration.

This module intentionally does not run Phaser yet.  It freezes and validates
the inputs that the Phaser execution layer will consume next.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Collection, Mapping

from .automr_input import (
    AutoMRInputError,
    ResolvedAutoMRInput,
    format_intent,
    read_intent,
    resolve_automr_input,
)
from .model_compatibility import (
    build_model_compatibility_facts,
    freeze_model_compatibility_facts,
)
from .model_assessment import (
    ModelAssessment,
    ModelAssessmentError,
    copy_preserving_model,
    file_sha256,
    inspect_pdb,
)
from .run_context import artifact_reference
from .phosphate import phosphate_intent_summary, validate_phosphate_intent, PhosphateError
from .backbone import make_backbone_policy
from .sequence_reference import SequenceReferenceError
from .sequence_family import prepare_sequence_family, freeze_sequence_family
from .search_model_comparison import freeze_search_model_comparison
from .frame_postmr import frame_postmr_spec
from .symmetry import StandardSymmetryAssessment, SymmetryError, assess_standard_symmetry
from .reflection_data import ReflectionDataError, describe, resolve_phaser_data


@dataclass(frozen=True)
class AutoMRPreflightResult:
    status: str
    message: str
    run_directory: Path
    generated_config: bool
    report_path: Path
    phosphate_summary: str | None = None


def _next_run_directory(dataset: Path) -> Path:
    parent = dataset / "AutoMR"
    parent.mkdir(exist_ok=True)
    for number in range(1, 10000):
        candidate = parent / f"run_{number:03d}"
        try:
            candidate.mkdir()
        except FileExistsError:
            continue
        return candidate
    raise AutoMRInputError(f"Could not allocate another run directory beneath {parent}")


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _post_mr_plan(resolved: ResolvedAutoMRInput) -> dict[str, object]:
    pair = None
    if resolved.pair:
        model_pair = None
        if resolved.model_pair:
            model_pair = {
                "tokens": [resolved.model_pair[0].token, resolved.model_pair[1].token],
                "ligand_codes": [
                    resolved.model_pair[0].ligand_code,
                    resolved.model_pair[1].ligand_code,
                ],
            }
        pair = {
            "requested": resolved.pair_text,
            "ligand_codes": [resolved.pair[0].ligand_code, resolved.pair[1].ligand_code],
            "model_pair": model_pair,
            "exact_model_match": resolved.exact_pair_model,
            "mutation_required": not bool(resolved.exact_pair_model),
            "ordered_roles": ["first standard site", "second standard site"],
            "site_assignment": "pending standard-frame site manifest",
        }
    return {
        "allow_op3_sites": list(resolved.allow_op3_sites),
        "backbone_policy": make_backbone_policy(
            resolved.backbone_sites,
            allow_unreviewed=resolved.allow_unreviewed_backbone,
        ),
        **({"phosphate_intent": resolved.phosphate_intent} if resolved.phosphate_intent is not None else {}),
        "application_order": ["sequences", "standard_pair", "explicit_mutations"],
        **(
            {"sequence_thread": dict(resolved.sequence_thread)}
            if resolved.sequence_thread is not None else {}
        ),
        "sequences": dict(resolved.sequences),
        "standard_pair": pair,
        "mutations": {
            site: {"requested": ligand.token, "ligand_code": ligand.ligand_code}
            for site, ligand in resolved.mutations.items()
        },
        "execution": "not performed by AutoMR preflight",
    }


def _validate_edit_targets(
    resolved: ResolvedAutoMRInput,
    assessment: ModelAssessment,
) -> None:
    """Require sequence lengths and explicit mutation sites to match the model."""
    thread_sequences = (
        resolved.sequence_thread.get("sequences", {})
        if resolved.sequence_thread is not None else {}
    )
    for source, sequences in (("Thread sequence", thread_sequences), ("Sequence", resolved.sequences)):
        if not isinstance(sequences, Mapping):
            raise AutoMRInputError(f"{source} declarations are malformed")
        for chain, sequence in sequences.items():
            model_length = assessment.polymer_residues_by_chain.get(chain)
            if model_length is None:
                available = ", ".join(assessment.chains)
                raise AutoMRInputError(
                    f"{source} chain {chain!r} is absent from the MR model; available chains: {available}"
                )
            if not isinstance(sequence, str) or len(sequence) != model_length:
                raise AutoMRInputError(
                    f"{source} for chain {chain} has length "
                    f"{len(sequence) if isinstance(sequence, str) else 'invalid'}, "
                    f"but the MR model contains {model_length} polymer residues"
                )
    thread_sites = (
        resolved.sequence_thread.get("site_codes", {})
        if resolved.sequence_thread is not None else {}
    )
    if not isinstance(thread_sites, Mapping):
        raise AutoMRInputError("Thread site-chemistry declarations are malformed")
    for site in thread_sites:
        chain, residue = (part.strip() for part in site.split(":", 1))
        if residue not in assessment.polymer_residue_ids_by_chain.get(chain, []):
            raise AutoMRInputError(f"Thread chemistry target {site} does not exist in the MR model")
    for site in resolved.mutations:
        chain, residue = (part.strip() for part in site.split(":", 1))
        if residue not in assessment.polymer_residue_ids_by_chain.get(chain, []):
            raise AutoMRInputError(f"Mutation target {site} does not exist in the MR model")

    if (
        resolved.frame is not None
        and isinstance(resolved.model_provider, Mapping)
        and resolved.model_provider.get("kind") == "explicit-standard-model"
    ):
        try:
            standard_sites = frame_postmr_spec(resolved.frame.name).sites
        except KeyError:
            standard_sites = ()
        for site in standard_sites:
            if site.resid not in assessment.polymer_residue_ids_by_chain.get(site.chain, []):
                raise AutoMRInputError(
                    f"Explicit standard model is incompatible with frame {resolved.frame.name}: "
                    f"required PostMR site {site.text} is absent"
                )

    for site in resolved.allow_op3_sites:
        chain, residue = site.split(":", 1)
        if residue not in assessment.polymer_residue_ids_by_chain.get(chain, []):
            raise AutoMRInputError(f"5'-phosphate request {site} does not exist in the MR model")
    for site in resolved.backbone_sites:
        chain, residue = site.split(":", 1)
        if residue not in assessment.polymer_residue_ids_by_chain.get(chain, []):
            raise AutoMRInputError(f"Backbone chemistry site {site} does not exist in the MR model")


def _default_mirror_transformer(source: Path, destination: Path) -> Path:
    try:
        from restraints.mirror_pdb import mirror_pdb
    except ImportError as exc:
        raise AutoMRInputError(
            "--mirror requires a NARestraints installation with restraints.mirror_pdb"
        ) from exc
    try:
        output = Path(mirror_pdb(source, destination)).expanduser().resolve()
    except Exception as exc:
        raise AutoMRInputError(f"NARestraints could not mirror the MR model: {exc}") from exc
    if output != destination.resolve() or not output.is_file():
        raise AutoMRInputError(
            f"NARestraints did not create the expected mirrored model: {destination}"
        )
    return output


def _validate_mirror_inventory(
    source: ModelAssessment,
    mirrored: ModelAssessment,
) -> None:
    comparisons = (
        ("atom count", source.atom_count, mirrored.atom_count),
        ("polymer residue count", source.polymer_residue_count, mirrored.polymer_residue_count),
        ("HETATM count", source.heteroatom_count, mirrored.heteroatom_count),
        ("chain inventory", source.chains, mirrored.chains),
        (
            "polymer residues by chain",
            source.polymer_residues_by_chain,
            mirrored.polymer_residues_by_chain,
        ),
        (
            "polymer residue IDs by chain",
            source.polymer_residue_ids_by_chain,
            mirrored.polymer_residue_ids_by_chain,
        ),
    )
    changed = [
        f"{label}: {before!r} -> {after!r}"
        for label, before, after in comparisons
        if before != after
    ]
    if changed:
        raise AutoMRInputError(
            "Mirroring changed the model inventory unexpectedly: " + "; ".join(changed)
        )


def _log_text(
    resolved: ResolvedAutoMRInput,
    assessment: ModelAssessment,
    symmetry: StandardSymmetryAssessment | None,
    status: str,
    config_path: Path,
    generated: bool,
    phaser_data: Mapping[str, object] | None = None,
) -> str:
    pair_line = "none"
    if resolved.pair:
        pair_line = (
            f"{resolved.pair_text} -> "
            f"{resolved.pair[0].ligand_code}:{resolved.pair[1].ligand_code}"
        )
    model_pair_line = "none"
    if resolved.model_pair:
        model_pair_line = (
            f"{resolved.model_pair[0].token}:{resolved.model_pair[1].token} -> "
            f"{resolved.model_pair[0].ligand_code}:{resolved.model_pair[1].ligand_code}"
        )
    mutation_lines = [
        f"  {site}: {ligand.token} -> {ligand.ligand_code}"
        for site, ligand in resolved.mutations.items()
    ] or ["  none"]
    warning_lines = [f"  {warning}" for warning in assessment.warnings] or ["  none"]
    catalogue_warning_lines = [
        f"  {warning}" for warning in resolved.catalogue_warnings
    ] or ["  none"]
    space_group = symmetry.evidence.mtz_symbol if symmetry else "not gated (nonstandard mode)"
    copies = str(symmetry.mr_copies) if symmetry else "not set by standard-frame rule"
    symmetry_warning = symmetry.red_flag if symmetry and symmetry.red_flag else "none"
    data_line = (
        f"Phaser data: {phaser_data['requested']} - "
        f"{describe(str(phaser_data['requested']))} ({phaser_data['labin']})"
        if phaser_data else "Phaser data: not recorded"
    )
    return "\n".join([
        "NASolve AutoMR preflight",
        f"Status: {status}",
        "Phaser executed: no",
        f"Dataset: {resolved.dataset.root}",
        f"Input configuration: {config_path} ({'generated' if generated else 'existing'})",
        f"Mode: {resolved.mode}",
        f"Frame: {resolved.frame.name if resolved.frame else 'none'}",
        phosphate_intent_summary(resolved.allow_op3_sites, resolved.phosphate_intent),
        ("Backbone chemistry: standard phosphodiester" if not resolved.backbone_sites else
         "Backbone chemistry exceptions: " + ", ".join(
             f"{site}={mode}" for site, mode in resolved.backbone_sites.items()
         )),
        f"Authoritative space group: {space_group}",
        "Space-group sources compared: "
        + (", ".join(symmetry.evidence.compared_sources) if symmetry else "none (nonstandard mode)"),
        f"Data_1*.cif: {resolved.dataset.metadata.name if resolved.dataset.metadata else 'absent'}",
        f"summary.html: {'present' if resolved.dataset.summary else 'absent'}",
        f"Planned MR copies: {copies}",
        f"Symmetry red flag: {symmetry_warning}",
        data_line,
        f"Pair: {pair_line}",
        f"Model pair: {model_pair_line}",
        f"Exact pair model: {resolved.exact_pair_model}",
        f"MR model: {resolved.model}",
        f"Model source: {resolved.model_source}",
        f"Mirrored before MR: {'yes' if resolved.mirror else 'no'}",
        "Catalogue warnings:",
        *catalogue_warning_lines,
        f"Model SHA-256: {assessment.sha256}",
        "Polymer residues (including modified HETATM nucleotides): "
        f"{assessment.polymer_residue_count}",
        f"Modified polymer residues from HETATM: {assessment.modified_polymer_residue_count}",
        f"HETATM atom records preserved: {assessment.heteroatom_count}",
        "Explicit mutations:",
        *mutation_lines,
        "Assessment warnings:",
        *warning_lines,
        "",
        "Next layer: use this frozen run input to prepare and execute Phenix Phaser.",
        "",
    ])


def prepare_automr(
    dataset: Path,
    config_path: Path | None = None,
    frame_override: str | None = None,
    pair_override: str | None = None,
    frames_dir: Path | None = None,
    allow_p1_standard: bool = False,
    mirror: bool = False,
    model_override: str | None = None,
    mtz_dump_executable: Path | None = None,
    phenix_environment: Mapping[str, str] | None = None,
    environ: Mapping[str, str] | None = None,
    valid_ligand_codes: Collection[str] | None = None,
    mirror_transformer: Callable[[Path, Path], Path] | None = None,
    *,
    resolved_input: ResolvedAutoMRInput | None = None,
    on_run_allocated: Callable[[Path], None] | None = None,
    phaser_data: str | None = None,
) -> AutoMRPreflightResult:
    """Validate and freeze one AutoMR request without executing Phaser.

    Campaigns may supply their checksum-verified, frozen selection. The same
    symmetry, model, and edit-target checks still run, but catalogue discovery
    and generation of a dataset input file do not run a second time.
    """
    root = dataset.expanduser().resolve()
    if resolved_input is None:
        selected_config = config_path.expanduser().resolve() if config_path else root / "nasolve.txt"
        existing_config = selected_config if selected_config.is_file() else None
        intent = read_intent(existing_config)
        resolved = resolve_automr_input(
            root,
            intent,
            frame_override=frame_override,
            pair_override=pair_override,
            frames_dir=frames_dir,
            allow_p1_standard=allow_p1_standard,
            mirror_override=mirror,
            model_override=model_override,
            environ=environ,
            valid_ligand_codes=valid_ligand_codes,
            phaser_data_override=phaser_data,
        )
    else:
        overrides = (
            config_path, frame_override, pair_override, frames_dir, model_override,
            phaser_data,
        )
        if any(value is not None for value in overrides) or allow_p1_standard or mirror:
            raise AutoMRInputError("A frozen AutoMR selection cannot be combined with input overrides")
        if resolved_input.dataset.root.resolve() != root:
            raise AutoMRInputError("Frozen AutoMR selection belongs to a different dataset")
        resolved = resolved_input
        selected_config = resolved.config_source or root / "nasolve.txt"
        existing_config = resolved.config_source
    symmetry: StandardSymmetryAssessment | None = None
    if resolved.mode == "standard":
        if mtz_dump_executable is None:
            raise AutoMRInputError(
                "Standard mode requires phenix.mtz.dump for the authoritative MTZ symmetry check"
            )
        try:
            symmetry = assess_standard_symmetry(
                resolved.dataset.reflections,
                resolved.dataset.metadata,
                resolved.dataset.summary,
                mtz_dump_executable,
                environment=phenix_environment,
                allow_p1_standard=resolved.allow_p1_standard,
            )
        except SymmetryError as exc:
            raise AutoMRInputError(str(exc)) from exc
    if resolved.phosphate_intent is not None:
        try:
            validate_phosphate_intent(resolved.phosphate_intent, resolved.allow_op3_sites)
        except PhosphateError as exc:
            raise AutoMRInputError(str(exc)) from exc
    source_assessment = inspect_pdb(
        resolved.model, polymer_ligand_codes=valid_ligand_codes
    )
    _validate_edit_targets(resolved, source_assessment)
    post_mr_plan = _post_mr_plan(resolved)
    family_seed = None
    if resolved.sequence_reference is not None:
        try:
            family_seed = prepare_sequence_family(
                resolved.sequence_reference, resolved.model, source_assessment,
                post_mr_plan,
                frame=resolved.frame.name if resolved.frame is not None else None,
                mirror=resolved.mirror,
            )
        except SequenceReferenceError as exc:
            raise AutoMRInputError(f"Sequence-family target preparation failed: {exc}") from exc
    effective_text = format_intent(resolved)
    # After the model checks, before a run is allocated: the MTZ must have the
    # requested merged F or I columns.
    try:
        phaser_data_record = resolve_phaser_data(
            resolved.phaser_data,
            resolved.dataset.reflections,
            mtz_dump_executable,
            environment=phenix_environment,
        )
    except ReflectionDataError as exc:
        raise AutoMRInputError(str(exc)) from exc

    generated = existing_config is None and resolved_input is None
    if generated:
        if config_path is not None and selected_config.parent != root:
            raise AutoMRInputError(
                "A missing --config path is not created outside the dataset; omit --config to generate nasolve.txt"
            )
        selected_config.write_text(effective_text, encoding="utf-8")

    run_dir = _next_run_directory(root)
    if on_run_allocated is not None:
        on_run_allocated(run_dir)
    model_dir = run_dir / "Model"
    model_dir.mkdir()
    copied_model = model_dir / "input_model.pdb"
    if resolved.mirror:
        source_copy = model_dir / "source_model.pdb"
        copy_preserving_model(resolved.model, source_copy, source_assessment)
        (mirror_transformer or _default_mirror_transformer)(source_copy, copied_model)
        assessment = inspect_pdb(
            copied_model, polymer_ligand_codes=valid_ligand_codes
        )
        _validate_mirror_inventory(source_assessment, assessment)
        assessment.heteroatoms_preserved = True
        assessment.copied_model = str(copied_model.resolve())
        _validate_edit_targets(resolved, assessment)
    else:
        assessment = source_assessment
        copy_preserving_model(resolved.model, copied_model, assessment)
    frame_sequence: Path | None = None
    if resolved.frame is not None:
        source_sequence = resolved.frame_sequence_source
        if source_sequence is None:
            legacy_sequence = resolved.model.parent / "seq_base.txt"
            source_sequence = legacy_sequence if legacy_sequence.is_file() else None
        if source_sequence is not None and source_sequence.is_file():
            frame_sequence = model_dir / "seq_base.txt"
            shutil.copyfile(source_sequence, frame_sequence)
    search_model_comparison = None
    comparison_payload = None
    if family_seed is not None:
        post_mr_plan["sequence_family"] = freeze_sequence_family(family_seed, run_dir)
        comparison_payload = json.loads(family_seed.comparison_bytes)
        search_model_comparison = freeze_search_model_comparison(
            comparison_payload, run_dir
        )
        post_mr_plan["application_order"] = [
            "sequence_family_reference", "thread_sequences", "sequences",
            "thread_site_chemistry", "standard_pair", "explicit_mutations",
        ]
    compatibility_facts = build_model_compatibility_facts(
        source_assessment=source_assessment,
        effective_assessment=assessment,
        model_provider=resolved.model_provider,
        candidate_construct_family=resolved.model_family,
        mode=resolved.mode,
        recipient_frame=resolved.frame.name if resolved.frame is not None else None,
        mirror_transform_applied=resolved.mirror,
        comparison=comparison_payload,
        comparison_artifact=search_model_comparison,
        terminal_phosphate_sites=tuple(resolved.allow_op3_sites),
        phosphate_intent=resolved.phosphate_intent,
        backbone_policy=post_mr_plan["backbone_policy"],
        symmetry=symmetry.to_dict() if symmetry is not None else None,
    )
    frozen_compatibility_facts = freeze_model_compatibility_facts(
        compatibility_facts, run_dir
    )
    _write_json(model_dir / "assessment.json", assessment.to_dict())
    (run_dir / "nasolve.input.txt").write_text(effective_text, encoding="utf-8")

    pair_needs_edit = bool(resolved.pair and not resolved.exact_pair_model)
    needs_edit = bool(
        pair_needs_edit or resolved.sequences or resolved.mutations
        or (family_seed is not None and family_seed.differences)
    )
    if symmetry and symmetry.red_flag:
        status = "READY_WITH_RED_FLAG"
        message = symmetry.red_flag
    elif needs_edit:
        status = "READY_POST_MR_MUTATION"
        message = "AutoMR inputs are frozen; a post-MR sequence/mutation plan is recorded"
    else:
        status = "READY"
        message = "AutoMR inputs are frozen and ready for the Phaser execution layer"
    report = {
        "schema_version": 1,
        "workflow": "automr",
        "stage": "preflight",
        "status": status,
        "message": message,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": str(root),
        "configuration": {
            "source": str(selected_config),
            "generated": generated,
            "effective_snapshot": str(run_dir / "nasolve.input.txt"),
        },
        "inputs": {
            "reflections": str(resolved.dataset.reflections),
            "reflections_sha256": file_sha256(resolved.dataset.reflections),
            # Optional autoPROC files; None when the dataset does not have them.
            "metadata": (
                str(resolved.dataset.metadata)
                if resolved.dataset.metadata is not None else None
            ),
            "summary": (
                str(resolved.dataset.summary)
                if resolved.dataset.summary is not None else None
            ),
            "model": str(resolved.model),
            "model_sha256": source_assessment.sha256,
            "model_selector": resolved.model_selector,
            "model_provider": resolved.model_provider,
            "mirror": resolved.mirror,
            "mirror_source_model": (
                str(model_dir / "source_model.pdb") if resolved.mirror else None
            ),
            "mirrored_model": str(copied_model) if resolved.mirror else None,
            "sequence_file": (
                str(resolved.sequence_file) if resolved.sequence_file else None
            ),
            "sequence_file_sha256": (
                file_sha256(resolved.sequence_file) if resolved.sequence_file else None
            ),
            **({"sequence_reference": post_mr_plan["sequence_family"]["reference"]}
               if family_seed is not None else {}),
            "frame_sequence": (
                {
                    **artifact_reference(frame_sequence, run_dir),
                    "sha256": file_sha256(frame_sequence),
                }
                if frame_sequence is not None
                else None
            ),
            "model_source": resolved.model_source,
            "catalogue_warnings": list(resolved.catalogue_warnings),
        },
        "mode": resolved.mode,
        "frame": asdict(resolved.frame) if resolved.frame else None,
        "phaser_data": phaser_data_record,
        "symmetry": symmetry.to_dict() if symmetry else {
            "gate": "not applied in nonstandard mode"
        },
        "model_assessment": assessment.to_dict(),
        "model_compatibility_facts": frozen_compatibility_facts,
        "post_mr_plan": post_mr_plan,
        **(
            {"search_model_comparison": search_model_comparison}
            if search_model_comparison is not None else {}
        ),
        **({
            "sequence_family_preflight": {
                "target_count": len(json.loads(family_seed.target_bytes)["sites"]),
                "scope": "literal-source-inventory-before-any-mirroring",
                "comparison_status": comparison_payload["status"],
                "differences": [
                    {"site": site, "before": before, "after": after}
                    for site, before, after in family_seed.differences
                ],
            },
        } if family_seed is not None else {}),
        "execution": {"phaser_ran": False},
    }
    report_path = run_dir / "report.json"
    _write_json(report_path, report)
    (run_dir / "automr.log").write_text(
        _log_text(
            resolved, assessment, symmetry, status, selected_config, generated,
            phaser_data_record,
        ),
        encoding="utf-8",
    )
    return AutoMRPreflightResult(
        status, message, run_dir, generated, report_path,
        phosphate_intent_summary(resolved.allow_op3_sites, resolved.phosphate_intent),
    )


__all__ = [
    "AutoMRInputError",
    "ModelAssessmentError",
    "AutoMRPreflightResult",
    "prepare_automr",
]
