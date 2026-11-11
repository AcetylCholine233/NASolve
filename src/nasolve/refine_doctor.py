"""Bounded, checkpoint-preserving refinement diagnosis and triage."""

from __future__ import annotations

import json
import math
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Mapping, Sequence

from .autorefine import (
    AutoRefineError,
    AutoRefineResult,
    build_reflection_plan,
    calculated_anomalous_groups,
    execute_autorefine,
    reflection_selector_policy,
    validate_refined_model,
    with_postmr_metal_sites,
    write_terminal_phosphate_protection as _write_terminal_phosphate_protection,
)
from .checkpoints import (
    CheckpointError,
    inherited_paths,
    initialize_registry,
    resolve_checkpoint,
)


class RefineDoctorError(RuntimeError):
    """Raised when refinement triage cannot be performed safely."""


@dataclass(frozen=True)
class FreeRAudit:
    status: str
    valid: bool | None
    details: dict[str, object]
    warnings: tuple[str, ...]
    log_path: Path


@dataclass(frozen=True)
class RefineDoctorResult:
    status: str
    message: str
    exit_code: int
    run_directory: Path
    doctor_directory: Path
    source_checkpoint: str
    current_checkpoint_preserved: bool
    recommended_checkpoint: str | None
    recommendation: str
    audit: FreeRAudit
    trials: tuple[AutoRefineResult, ...]
    benchmark: tuple[dict[str, object], ...]
    report_path: Path


@dataclass(frozen=True)
class RefineDoctorTrial:
    """One declarative branch; project presets may supply additional recipes later."""

    recipe: str
    use_experimental_phases: bool
    real_space_sites: bool = False
    adp_mode: str = "group"
    refine_occupancies: bool = False
    anomalous_mode: str = "off"
    requires_phases: bool = False
    maximum_resolution_limit: float | None = None
    minimum_observations_per_atom: float | None = None
    refine_coordinates: bool = True


DEFAULT_TRIALS = (
    RefineDoctorTrial("RefineDoctor/ML-group", False),
    RefineDoctorTrial(
        "RefineDoctor/MLHL-group", True, requires_phases=True
    ),
    RefineDoctorTrial(
        "RefineDoctor/ML-individual-ADP",
        False,
        adp_mode="individual",
        maximum_resolution_limit=3.2,
        minimum_observations_per_atom=3.0,
    ),
    RefineDoctorTrial("RefineDoctor/ML-coordinates-only", False, adp_mode="none"),
    RefineDoctorTrial("RefineDoctor/ML-group-B-only", False, refine_coordinates=False),
)

ANOMALOUS_TRIALS = (
    RefineDoctorTrial("RefineDoctor/ML-fixed-scattering", False, anomalous_mode="fixed"),
    RefineDoctorTrial("RefineDoctor/MLHL-fixed-scattering", True,
                     anomalous_mode="fixed", requires_phases=True),
    RefineDoctorTrial("RefineDoctor/ML-fdp-only", False, anomalous_mode="fdp-only"),
    RefineDoctorTrial("RefineDoctor/ML-coordinates-only-fixed-scattering", False,
                     adp_mode="none", anomalous_mode="fixed"),
    RefineDoctorTrial("RefineDoctor/ML-group-B-only-fixed-scattering", False,
                     refine_coordinates=False, anomalous_mode="fixed"),
)


_AUDIT_MARKER = "NASOLVE_FREE_R_AUDIT_JSON:"
_AUDIT_SCRIPT = r'''
from __future__ import print_function
import json
import sys
from iotbx import reflection_file_reader

# NASOLVE_FREE_R_AUDIT
filename, wanted_label, test_text = sys.argv[1:4]
test_value = int(test_text)
reflection_file = reflection_file_reader.any_reflection_file(filename)
try:
    arrays = reflection_file.as_miller_arrays(merge_equivalents=False)
except TypeError:
    arrays = reflection_file.as_miller_arrays()
selected = None
available = []
for array in arrays:
    info = array.info()
    labels = list(getattr(info, "labels", []) or [])
    label_string = info.label_string()
    available.append(label_string)
    if wanted_label in labels or label_string == wanted_label:
        selected = array
        break
if selected is None:
    raise RuntimeError("Free-R label not found; available: " + "; ".join(available))

unit_cell = selected.unit_cell()
indices = list(selected.indices())
data = list(selected.data())
groups = []
matching_method = "cctbx-match-bijvoet-mates"
try:
    matches = selected.match_bijvoet_mates()
    paired_positions = set()
    for first, second in matches.pairs():
        first, second = int(first), int(second)
        paired_positions.update([first, second])
        groups.append({
            "indices": [tuple(indices[first]), tuple(indices[second])],
            "values": set([int(data[first]), int(data[second])]),
            "d": float(unit_cell.d(indices[first])),
        })
    for position in range(len(indices)):
        if position in paired_positions:
            continue
        groups.append({
            "indices": [tuple(indices[position])],
            "values": set([int(data[position])]),
            "d": float(unit_cell.d(indices[position])),
        })
except Exception:
    matching_method = "exact-hkl-fallback"
    exact = {}
    for hkl, raw_value in zip(indices, data):
        index = tuple(int(value) for value in hkl)
        inverse = tuple(-value for value in index)
        key = min(index, inverse)
        item = exact.setdefault(key, {"indices": set(), "values": set(), "d": None})
        item["indices"].add(index)
        item["values"].add(int(raw_value))
        if item["d"] is None:
            item["d"] = float(unit_cell.d(index))
    groups = list(exact.values())

ordered = sorted(groups, key=lambda item: item["d"], reverse=True)
shells = [{"groups": 0, "free_groups": 0} for unused in range(10)]
free_groups = 0
inconsistent = 0
paired = 0
for number, item in enumerate(ordered):
    is_free = test_value in item["values"]
    if is_free:
        free_groups += 1
    if len(item["values"]) > 1:
        inconsistent += 1
    if len(item["indices"]) > 1:
        paired += 1
    shell = min(9, int(number * 10 / max(1, len(ordered))))
    shells[shell]["groups"] += 1
    shells[shell]["free_groups"] += int(is_free)

payload = {
    "array_labels": selected.info().label_string(),
    "array_anomalous": bool(selected.anomalous_flag()),
    "matching_method": matching_method,
    "stored_observations": len(selected.indices()),
    "independent_friedel_groups": len(ordered),
    "paired_friedel_groups": paired,
    "free_independent_groups": free_groups,
    "free_fraction": float(free_groups) / max(1, len(ordered)),
    "inconsistent_friedel_flag_groups": inconsistent,
    "test_flag_value": test_value,
    "resolution_shells": shells,
}
print("NASOLVE_FREE_R_AUDIT_JSON:" + json.dumps(payload, sort_keys=True))
'''


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _next_doctor_directory(run: Path) -> Path:
    root = run / "RefineDoctor"
    root.mkdir(parents=True, exist_ok=True)
    numbers = [
        int(match.group(1))
        for path in root.glob("doctor_*")
        if path.is_dir() and (match := re.fullmatch(r"doctor_(\d+)", path.name))
    ]
    destination = root / f"doctor_{max(numbers, default=0) + 1:03d}"
    destination.mkdir()
    return destination


def audit_free_r_flags(
    observations: Path,
    free_r_label: str,
    test_flag_value: int,
    phenix_python: Path,
    destination: Path,
    environment: Mapping[str, str] | None = None,
) -> FreeRAudit:
    """Audit the frozen test set without changing or regenerating it."""
    script_path = destination / "audit_free_r.py"
    log_path = destination / "free_r_audit.log"
    script_path.write_text(_AUDIT_SCRIPT.lstrip(), encoding="utf-8")
    try:
        completed = subprocess.run(
            [
                str(phenix_python.expanduser().resolve()),
                str(script_path),
                str(observations),
                free_r_label,
                str(test_flag_value),
            ],
            env=dict(environment) if environment is not None else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=180,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        log_path.write_text(f"Free-R audit could not run: {exc}\n", encoding="utf-8")
        return FreeRAudit(
            "UNAVAILABLE", None, {}, ("Free-R audit could not be executed",), log_path
        )
    log_path.write_text(completed.stdout, encoding="utf-8")
    matches = re.findall(
        re.escape(_AUDIT_MARKER) + r"(\{.*\})", completed.stdout, re.MULTILINE
    )
    if completed.returncode or not matches:
        return FreeRAudit(
            "UNAVAILABLE",
            None,
            {"phenix_python_exit_code": completed.returncode},
            ("Free-R audit returned no machine-readable result",),
            log_path,
        )
    try:
        details = json.loads(matches[-1])
    except json.JSONDecodeError:
        return FreeRAudit(
            "UNAVAILABLE", None, {}, ("Free-R audit result was malformed",), log_path
        )
    if not isinstance(details, dict):
        return FreeRAudit("UNAVAILABLE", None, {}, ("Free-R audit result was not an object",), log_path)

    warnings: list[str] = []
    invalid: list[str] = []
    inconsistent = details.get("inconsistent_friedel_flag_groups")
    if isinstance(inconsistent, int) and inconsistent:
        invalid.append(f"{inconsistent} Friedel group(s) have inconsistent test flags")
    fraction = details.get("free_fraction")
    if not isinstance(fraction, (int, float)) or not 0.02 <= float(fraction) <= 0.20:
        invalid.append("Free-R fraction is outside the guarded 2-20% range")
    free_groups = details.get("free_independent_groups")
    if not isinstance(free_groups, int) or free_groups < 10:
        invalid.append("Fewer than 10 independent Free-R groups are available")
    elif free_groups < 100:
        warnings.append(
            f"Only {free_groups} independent Free-R groups are available; Rfree is noisy"
        )
    paired = details.get("paired_friedel_groups")
    if paired == 0:
        warnings.append("Stored array did not expose Friedel mates for a direct consistency test")
    shells = details.get("resolution_shells")
    if isinstance(shells, list) and any(
        isinstance(shell, Mapping) and shell.get("groups", 0) and shell.get("free_groups") == 0
        for shell in shells
    ):
        warnings.append("At least one resolution shell contains no independent Free-R group")
    if invalid:
        return FreeRAudit("INVALID", False, details, tuple(invalid + warnings), log_path)
    if warnings:
        return FreeRAudit("NOISY", True, details, tuple(warnings), log_path)
    return FreeRAudit("VALID", True, details, (), log_path)


def _coordinate_atom_count(model: Path) -> int:
    return sum(
        line[0:6].strip().upper() in {"ATOM", "HETATM"}
        for line in model.read_text(encoding="utf-8", errors="replace").splitlines()
    )


def _benchmark_from_checkpoint(
    checkpoint: Mapping[str, object], source: str
) -> list[dict[str, object]]:
    metrics = checkpoint.get("metrics")
    values = metrics.get("anomalous_scatterers") if isinstance(metrics, Mapping) else None
    if not isinstance(values, list):
        return []
    benchmarks: list[dict[str, object]] = []
    for value in values:
        if not isinstance(value, Mapping) or not isinstance(
            value.get("refined_f_double_prime"), (int, float)
        ):
            continue
        benchmarks.append({**dict(value), "source_checkpoint": source})
    return benchmarks


def _metrics(checkpoint: Mapping[str, object]) -> dict[str, object]:
    values = checkpoint.get("metrics")
    return dict(values) if isinstance(values, Mapping) else {}


def _terminal_geometry_audit(
    checkpoint: Mapping[str, object],
) -> dict[str, object] | None:
    metrics = checkpoint.get("metrics")
    audit = (
        metrics.get("terminal_geometry_audit")
        if isinstance(metrics, Mapping)
        else None
    )
    if isinstance(audit, Mapping) and audit.get("requires_review") is True:
        return dict(audit)
    return None


def _terminal_geometry_rescued(
    source_audit: Mapping[str, object],
    trial_audit: object,
) -> bool:
    if not isinstance(trial_audit, Mapping):
        return False
    if (
        trial_audit.get("status") != "PASS"
        or trial_audit.get("requires_review") is not False
    ):
        return False

    source_rows = source_audit.get("sites")
    trial_rows = trial_audit.get("sites")
    if not isinstance(source_rows, list) or not isinstance(trial_rows, list):
        return False

    source_sites = {
        row.get("site")
        for row in source_rows
        if isinstance(row, Mapping)
        and row.get("requires_review") is True
        and isinstance(row.get("site"), str)
    }
    trial_sites = {
        row.get("site")
        for row in trial_rows
        if isinstance(row, Mapping)
        and row.get("requires_review") is False
        and isinstance(row.get("site"), str)
    }
    return bool(source_sites) and source_sites.issubset(trial_sites)


def write_terminal_phosphate_protection(
    audit: Mapping[str, object],
    destination: Path,
    *,
    sigma: float = 1.0,
) -> dict[str, object]:
    """Doctor-facing wrapper around AutoRefine terminal-geometry protection."""
    try:
        return _write_terminal_phosphate_protection(
            audit,
            destination,
            sigma=sigma,
        )
    except AutoRefineError as exc:
        raise RefineDoctorError(str(exc)) from exc

def _candidate(
    checkpoint: str, recipe: str, statistics: Mapping[str, object], *,
    usable: bool = True, accepted: bool = True,
) -> dict[str, object]:
    def metric(key: str) -> float | None:
        value = statistics.get(key)
        return float(value) if (
            not isinstance(value, bool) and isinstance(value, (int, float))
            and math.isfinite(value) and 0 <= value <= 1
        ) else None
    work, free = metric("r_work"), metric("r_free")
    terminal_audit = statistics.get("terminal_geometry_audit")
    terminal_geometry_review = (
        isinstance(terminal_audit, Mapping)
        and terminal_audit.get("requires_review") is True
    )
    eligible = (
        usable
        and work is not None
        and free is not None
        and not terminal_geometry_review
    )
    gap = free - work if work is not None and free is not None else None
    return {
        "checkpoint": checkpoint,
        "recipe": recipe,
        "r_work": work,
        "r_free": free,
        "r_free_minus_r_work": gap,
        "usable": usable,
        "terminal_geometry_review": terminal_geometry_review,
        "eligible_for_inspection": eligible,
        "strict_success": bool(eligible and accepted and work < .30 and work < free),
    }


def _rank_for_inspection(candidates: list[dict[str, object]]) -> list[dict[str, object]]:
    """Descriptive numerical order, with no claim of statistical significance."""
    return sorted(
        (item for item in candidates if item.get("eligible_for_inspection")),
        key=lambda item: (item["r_free"], item["r_work"]),
    )


def _recommend(
    candidates: list[dict[str, object]],
    source_checkpoint: str,
    audit: FreeRAudit,
) -> tuple[str, str | None, str, int]:
    source = candidates[0]
    if audit.valid is not True:
        return "REFINE_DOCTOR_REVIEW", None, "Free-R integrity must be established before a recommendation", 2
    if source.get("strict_success") is True:
        return (
            "REFINE_DOCTOR_GOOD_ENOUGH",
            source_checkpoint,
            "The selected checkpoint already meets the numerical gate; no replacement is needed",
            0,
        )
    successes = [candidate for candidate in candidates[1:] if candidate.get("strict_success")]
    if successes:
        choice = successes[0]
        return (
            "REFINE_DOCTOR_RECOMMEND",
            str(choice["checkpoint"]),
            "The first eligible branch passed the numerical gate; map/model inspection is still required",
            0,
        )
    ranked = _rank_for_inspection(candidates)
    detail = (
        f" Lowest Rfree among usable candidates: {ranked[0]['checkpoint']}; "
        "this numerical ordering does not establish superiority given test-set uncertainty."
        if ranked else " No usable candidate has valid comparison statistics."
    )
    return (
        "REFINE_DOCTOR_REVIEW",
        None,
        "No branch passed the numerical gate; inspect the model and maps." + detail,
        2,
    )


def _update_run_report(run: Path, payload: Mapping[str, object]) -> None:
    path = run / "report.json"
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RefineDoctorError(f"Could not update run report: {exc}") from exc
    history = report.setdefault("refine_doctor_history", [])
    if not isinstance(history, list):
        raise RefineDoctorError("Run report has malformed Refine Doctor history")
    history.append(dict(payload))
    report["refine_doctor"] = dict(payload)
    report["stage"] = "refine-doctor"
    report["status"] = payload["status"]
    report["message"] = payload["message"]
    report["updated_utc"] = payload["created_utc"]
    _write_json(path, report)


def execute_refine_doctor(
    run_directory: Path,
    refine_executable: Path,
    mtz_dump_executable: Path,
    *,
    phenix_version: str,
    environment: Mapping[str, str] | None = None,
    from_checkpoint: str | None = None,
    macro_cycles: int = 3,
    processor_count: int | None = None,
    trial_recipes: Sequence[RefineDoctorTrial] | None = None,
    max_trials: int = 5,
    progress: Callable[[str, str, Path], None] | None = None,
) -> RefineDoctorResult:
    """Audit and run a finite set of sibling refinements without selecting one."""
    if isinstance(max_trials, bool) or not isinstance(max_trials, int) or not 1 <= max_trials <= 10:
        raise RefineDoctorError("Doctor max_trials must be an integer from 1 to 10")
    if isinstance(macro_cycles, bool) or not isinstance(macro_cycles, int) or not 1 <= macro_cycles <= 10:
        raise RefineDoctorError("Doctor macrocycles must be an integer from 1 to 10")
    try:
        selector_policy = reflection_selector_policy(phenix_version)
    except AutoRefineError as exc:
        raise RefineDoctorError(str(exc)) from exc
    try:
        run, registry = initialize_registry(run_directory)
        source = resolve_checkpoint(registry, from_checkpoint)
        if source.get("usable") is not True:
            raise RefineDoctorError(f"Checkpoint {source.get('id')} is not reusable")
        inherited = inherited_paths(source, run)
    except CheckpointError as exc:
        raise RefineDoctorError(str(exc)) from exc
    source_id = str(source["id"])
    original_current = str(registry.get("current"))
    source_terminal_audit = _terminal_geometry_audit(source)
    observations = inherited.get("observations")
    model = inherited.get("model")
    if not isinstance(observations, Path) or not isinstance(model, Path):
        raise RefineDoctorError("Selected checkpoint has incomplete crystallographic inputs")
    try:
        report = with_postmr_metal_sites(
            json.loads((run / "report.json").read_text(encoding="utf-8")), run
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise RefineDoctorError(f"Could not read run report: {exc}") from exc
    except AutoRefineError as exc:
        raise RefineDoctorError(str(exc)) from exc
    try:
        validate_refined_model(model, report)
    except AutoRefineError as exc:
        raise RefineDoctorError(str(exc)) from exc
    destination = _next_doctor_directory(run)
    try:
        plan = build_reflection_plan(
            report, observations, mtz_dump_executable, environment, phase_file=None
        )
    except AutoRefineError as exc:
        raise RefineDoctorError(str(exc)) from exc
    phenix_python = refine_executable.expanduser().resolve().parent / "phenix.python"
    audit = audit_free_r_flags(
        observations,
        plan.free_r_label,
        plan.free_r_test_value,
        phenix_python,
        destination,
        environment,
    )
    created = _now()
    if audit.valid is not True:
        status = "REFINE_DOCTOR_FLAG_REPAIR_REQUIRED" if audit.valid is False else "REFINE_DOCTOR_REVIEW"
        message = "Free-R audit failed or was unavailable; no refinement trials were run"
        recommendation = (
            "Repair the flags deliberately in a separate branch; NASolve did not regenerate them"
            if audit.valid is False else "Restore the Phenix Free-R audit before retrying Doctor"
        )
        payload = {
            "status": status,
            "message": message,
            "created_utc": created,
            "phenix_version": selector_policy.phenix_version,
            "reflection_selector_mode": selector_policy.mode,
            "source_checkpoint": source_id,
            "current_checkpoint_preserved": True,
            "audit": {**audit.details, "status": audit.status, "warnings": list(audit.warnings)},
            "benchmark": [],
            "trials": [],
            "recommended_checkpoint": None,
            "recommendation": recommendation,
            "inspection_checkpoint": None,
            "automatic_selection": False,
            "triage": {"stop_reason": "free-r-audit", "max_trials": max_trials,
                       "recipes": [], "next_actions": [recommendation]},
        }
        report_path = destination / "report.json"
        _write_json(report_path, payload)
        _update_run_report(run, payload)
        return RefineDoctorResult(
            status, message, 2, run, destination, source_id, True, None,
            str(payload["recommendation"]), audit, (), (), report_path
        )

    terminal_parent_id: str | None = None
    terminal_protection_path: Path | None = None
    terminal_protection: dict[str, object] | None = None

    source_metrics = _metrics(source)
    if source_terminal_audit is not None:
        source_metrics["terminal_geometry_audit"] = source_terminal_audit

    if source_terminal_audit is not None and trial_recipes is None:
        parent_value = source.get("parent")
        if not isinstance(parent_value, str):
            raise RefineDoctorError(
                "Terminal-geometry rescue requires a reusable parent checkpoint"
            )
        try:
            clean_parent = resolve_checkpoint(registry, parent_value)
        except CheckpointError as exc:
            raise RefineDoctorError(str(exc)) from exc
        if clean_parent.get("usable") is not True:
            raise RefineDoctorError(
                f"Terminal-geometry rescue parent {parent_value} is not reusable"
            )
        if _terminal_geometry_audit(clean_parent) is not None:
            raise RefineDoctorError(
                f"Terminal-geometry rescue parent {parent_value} is itself "
                "chemically flagged; refusing to call it a clean parent"
            )

        terminal_parent_id = parent_value
        terminal_protection_path = (
            destination / "terminal_phosphate_protection.phil"
        )
        terminal_protection = write_terminal_phosphate_protection(
            source_terminal_audit,
            terminal_protection_path,
        )

    source_candidate = _candidate(
        source_id,
        str(source.get("recipe", "source")),
        source_metrics,
        usable=source.get("usable") is True,
        accepted=source.get("status") == "SUCCESS",
    )
    candidates = [source_candidate]
    benchmark = _benchmark_from_checkpoint(source, source_id)
    trials: list[AutoRefineResult] = []
    phase_available = isinstance(inherited.get("phases"), Path)
    autosol = report.get("autosol")
    if isinstance(autosol, Mapping) and autosol.get("use_for_refinement") is True:
        phase_available = True
    independent = audit.details.get("independent_friedel_groups")
    ratio = float(independent) / max(1, _coordinate_atom_count(model)) if isinstance(independent, int) else None
    pool = (
        tuple(trial_recipes)
        if trial_recipes is not None
        else (
            RefineDoctorTrial(
                "RefineDoctor/terminal-phosphate-protected",
                True,
                real_space_sites=True,
                adp_mode="group",
                refine_occupancies=True,
                anomalous_mode="refine",
            ),
        )
        if source_terminal_audit is not None
        else ANOMALOUS_TRIALS
        if plan.anomalous
        else DEFAULT_TRIALS
    )
    stop_reason = (
        "source-passes" if source_candidate["strict_success"] and trial_recipes is None
        else "anomalous-data-unavailable" if plan.anomalous_fallback and phase_available
        else None
    )
    scattering = None
    scattering_error = None
    if stop_reason is None and plan.anomalous and any(
        spec.anomalous_mode in {"fixed", "fdp-only"} for spec in pool
    ):
        try:
            scattering = calculated_anomalous_groups(report, refine_executable, environment)
        except AutoRefineError as exc:
            scattering_error = str(exc)
    recipe_records: list[dict[str, object]] = []
    trial_specs: list[RefineDoctorTrial] = []
    for spec in pool:
        record = {"recipe": spec.recipe, "state": "skipped", "reason": ""}
        recipe_records.append(record)
        reason = None
        if stop_reason is not None:
            reason = stop_reason
        elif spec.requires_phases and not phase_available:
            reason = "Accepted experimental phases are unavailable"
        elif spec.maximum_resolution_limit is not None and (
            plan.resolution_limit is None or plan.resolution_limit > spec.maximum_resolution_limit
        ):
            reason = "The resolution does not support this ADP recipe"
        elif spec.minimum_observations_per_atom is not None and (
            ratio is None or ratio < spec.minimum_observations_per_atom
        ):
            reason = "Too few independent observations per atom"
        elif spec.anomalous_mode in {"fixed", "fdp-only"} and scattering is None:
            reason = scattering_error or "Explicit scattering requires anomalous data and wavelength"
        elif len(trials) >= max_trials:
            stop_reason = reason = "trial-budget-exhausted"
        if reason is not None:
            record["reason"] = reason
            continue
        record.update(state="running", reason="Eligible in declared recipe order")
        trial_specs.append(spec)
        recipe = spec.recipe

        def trial_progress(checkpoint: str, log: Path, recipe_name: str = recipe) -> None:
            if progress is not None:
                progress(recipe_name, checkpoint, log)

        terminal_trial = (
            recipe == "RefineDoctor/terminal-phosphate-protected"
        )
        trial_parent = (
            terminal_parent_id
            if terminal_trial
            else source_id
        )
        trial_extra_restraints = (
            (terminal_protection_path,)
            if terminal_trial and terminal_protection_path is not None
            else ()
        )

        if trial_parent is None:
            record.update(
                state="failed",
                reason="Terminal-geometry rescue parent was not resolved",
            )
            stop_reason = "technical-failure"
            continue

        try:
            result = execute_autorefine(
                run,
                refine_executable,
                mtz_dump_executable,
                phenix_version=selector_policy.phenix_version,
                environment=environment,
                from_checkpoint=trial_parent,
                recipe=recipe,
                macro_cycles=macro_cycles,
                processor_count=processor_count,
                use_experimental_phases=spec.use_experimental_phases,
                real_space_sites=spec.real_space_sites,
                adp_mode=spec.adp_mode,
                refine_occupancies=spec.refine_occupancies,
                anomalous_mode=spec.anomalous_mode,
                refine_coordinates=spec.refine_coordinates,
                anomalous_groups=scattering if spec.anomalous_mode in {"fixed", "fdp-only"} else None,
                extra_restraints=trial_extra_restraints,
                terminal_geometry_protection=(
                    terminal_protection if terminal_trial else None
                ),
                auto_select_success=False,
                progress=trial_progress,
            )
        except AutoRefineError as exc:
            record.update(state="failed", reason=str(exc))
            stop_reason = "technical-failure"
            continue
        trials.append(result)
        record.update(state="completed", checkpoint=result.checkpoint_id, reason=result.message)
        usable = result.status in {"AUTOREFINE_READY", "AUTOREFINE_REVIEW"} and result.model_path is not None
        candidates.append(_candidate(result.checkpoint_id, recipe, result.statistics,
                                      usable=usable, accepted=result.status == "AUTOREFINE_READY"))
        if usable and not benchmark and spec.anomalous_mode in {"refine", "fdp-only"}:
            benchmark = [
                {**dict(item), "source_checkpoint": result.checkpoint_id}
                for item in result.statistics.get("anomalous_scatterers", [])
                if isinstance(item, Mapping)
                and isinstance(item.get("refined_f_double_prime"), (int, float))
            ]
        if not usable:
            stop_reason = "technical-failure"
        elif terminal_trial:
            trial_terminal = result.statistics.get(
                "terminal_geometry_audit"
            )
            if (
                source_terminal_audit is not None
                and _terminal_geometry_rescued(
                    source_terminal_audit,
                    trial_terminal,
                )
            ):
                stop_reason = "terminal-geometry-rescued"
            else:
                stop_reason = "terminal-geometry-unresolved"
        elif candidates[-1]["strict_success"]:
            stop_reason = "numerical-pass"
    if stop_reason is None:
        stop_reason = "eligible-recipes-exhausted"

    try:
        _, final_registry = initialize_registry(run)
    except CheckpointError as exc:
        raise RefineDoctorError(str(exc)) from exc
    current_preserved = str(final_registry.get("current")) == original_current
    if not current_preserved:
        raise RefineDoctorError("Refine Doctor changed the current checkpoint unexpectedly")
    status, recommended, recommendation, exit_code = _recommend(
        candidates, source_id, audit
    )
    message = {
        "REFINE_DOCTOR_GOOD_ENOUGH": "The source meets the numerical gate; model/map inspection remains required",
        "REFINE_DOCTOR_RECOMMEND": "A bounded refinement branch is recommended for inspection",
        "REFINE_DOCTOR_REVIEW": "Refinement remains a user-review case after bounded triage",
    }[status]
    payload = {
        "status": status,
        "message": message,
        "created_utc": created,
        "phenix_version": selector_policy.phenix_version,
        "reflection_selector_mode": selector_policy.mode,
        "source_checkpoint": source_id,
        "current_checkpoint_preserved": current_preserved,
        "audit": {**audit.details, "status": audit.status, "warnings": list(audit.warnings)},
        "benchmark": benchmark,
        "terminal_geometry": {
            "triggered": source_terminal_audit is not None,
            "source_audit": source_terminal_audit,
            "clean_parent_checkpoint": terminal_parent_id,
            "protection": terminal_protection,
        },
        "eligibility": {
            "resolution_limit": plan.resolution_limit,
            "independent_observations_per_atom": ratio,
            "individual_adp_trial": any(
                spec.recipe.endswith("individual-ADP") for spec in trial_specs
            ),
        },
        "candidates": candidates,
        "inspection_checkpoint": (
            _rank_for_inspection(candidates)[0]["checkpoint"]
            if _rank_for_inspection(candidates) else None
        ),
        "ranking": {
            "basis": "Rfree then Rwork, among usable comparable candidates",
            "statistical_superiority_established": False,
            "checkpoints": [item["checkpoint"] for item in _rank_for_inspection(candidates)],
        },
        "triage": {
            "observation_mode": "anomalous" if plan.anomalous else "mean",
            "max_trials": max_trials,
            "macrocycles_per_trial": macro_cycles,
            "stop_reason": stop_reason,
            "recipes": recipe_records,
            "next_actions": (
                [
                    "Inspect the locally protected terminal phosphate and maps before selecting it",
                    "Confirm the terminal-geometry audit remains below the review threshold",
                ]
                if stop_reason == "terminal-geometry-rescued"
                else ["Inspect the failed trial log or missing prerequisites before retrying"]
                if stop_reason in {
                    "technical-failure",
                    "anomalous-data-unavailable",
                    "terminal-geometry-unresolved",
                } or scattering_error
                else ["Inspect candidate maps, modified-site restraints, and B factors",
                      "Review completeness, anisotropy, scaling, and twinning evidence if fit remains poor",
                      "Keep the frozen Free-R flags; no statistically justified winner has been established"]
                if recommended is None
                else ["Inspect the recommended model and maps before selecting it"]
            ),
        },
        "trials": [
            {
                "checkpoint": trial.checkpoint_id,
                "parent": trial.parent_checkpoint,
                "status": trial.status,
                "report": str(trial.report_path),
                "log": str(trial.log_path),
            }
            for trial in trials
        ],
        "recommended_checkpoint": recommended,
        "recommendation": recommendation,
        "automatic_selection": False,
    }
    report_path = destination / "report.json"
    _write_json(report_path, payload)
    _update_run_report(run, payload)
    return RefineDoctorResult(
        status=status,
        message=message,
        exit_code=exit_code,
        run_directory=run,
        doctor_directory=destination,
        source_checkpoint=source_id,
        current_checkpoint_preserved=current_preserved,
        recommended_checkpoint=recommended,
        recommendation=recommendation,
        audit=audit,
        trials=tuple(trials),
        benchmark=tuple(benchmark),
        report_path=report_path,
    )


__all__ = [
    "FreeRAudit",
    "RefineDoctorError",
    "RefineDoctorResult",
    "RefineDoctorTrial",
    "DEFAULT_TRIALS",
    "audit_free_r_flags",
    "write_terminal_phosphate_protection",
    "execute_refine_doctor",
]
