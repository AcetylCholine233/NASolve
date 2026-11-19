"""Phenix Phaser execution and TFZ gating for a frozen AutoMR run."""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

from .model_assessment import file_sha256
from .reflection_data import (
    AMPLITUDES,
    INTENSITIES,
    describe,
    reported_data_kind,
    reported_labin,
)
from .run_context import resolve_artifact_path


class PhaserExecutionError(RuntimeError):
    """Raised when a frozen run cannot be prepared for Phaser."""


@dataclass(frozen=True)
class PhaserExecutionResult:
    status: str
    message: str
    run_directory: Path
    phaser_directory: Path
    report_path: Path
    log_path: Path
    parameter_path: Path
    tfz: float | None   # the TFZ that decided the status (after refinement when reported)
    llg: float | None
    solution_pdb: Path | None
    solution_mtz: Path | None
    translation_tfz: float | None = None   # TFZ= from the translation search
    refined_tfz: float | None = None       # TFZ== after refinement
    tfz_basis: str | None = None           # "refined" or "translation"

    @property
    def tfz_scores(self) -> "TFZScores":
        return TFZScores(self.tfz, self.llg, self.translation_tfz, self.refined_tfz, self.tfz_basis)

    @property
    def exit_code(self) -> int:
        if self.status == "MR_SUCCESS":
            return 0
        if self.status == "MR_REVIEW":
            return 3
        return 4


#: Best TFZ below this is MR_FAILED; at or above TFZ_THRESHOLD it is MR_SUCCESS;
#: between the two it is MR_REVIEW.
TFZ_FLOOR = 7.0
TFZ_THRESHOLD = 8.0


class TFZOverrideError(ValueError):
    """A run that ``--tfz-override`` cannot accept."""


def check_tfz_override(report: Mapping[str, object], run_directory: Path) -> dict[str, object] | None:
    """Accept a Phaser result whose only problem is a low TFZ (``--tfz-override``).

    The cut-offs do not change: the run keeps its MR_REVIEW or MR_FAILED
    status, and the acceptance is recorded by the step that uses it. A run
    qualifies when Phaser exited with status 0, reported a TFZ, wrote both
    the solution model and map, and used the reflection data that were
    requested. Returns None when the TFZ is at or above TFZ_THRESHOLD
    (nothing to accept), the record to store when it is below, and raises
    TFZOverrideError when the run failed for another reason.
    """
    execution = report.get("execution")
    phaser = execution.get("phaser") if isinstance(execution, Mapping) else None
    if not isinstance(phaser, Mapping):
        raise TFZOverrideError("The run has no Phaser record; --tfz-override needs a Phaser result")
    return_code = phaser.get("return_code")
    if return_code != 0:
        raise TFZOverrideError(
            f"Phaser exited with status {return_code}; --tfz-override only accepts a low TFZ"
        )
    tfz = phaser.get("best_tfz")
    if not isinstance(tfz, (int, float)) or isinstance(tfz, bool) or not math.isfinite(tfz):
        raise TFZOverrideError("Phaser reported no TFZ; --tfz-override only accepts a low TFZ")
    reflection = phaser.get("reflection_data")
    if isinstance(reflection, Mapping):
        requested, used = reflection.get("requested"), reflection.get("phaser_reported_kind")
        if requested is not None and used is not None and used != requested:
            raise TFZOverrideError(
                f"Phaser used {used} data although {requested} data were requested; "
                "--tfz-override only accepts a low TFZ"
            )
    for key, suffix, name in (("solution_pdb", ".pdb", "model"), ("solution_mtz", ".mtz", "map")):
        value = phaser.get(key)
        path = resolve_artifact_path(value, run_directory) if isinstance(value, str) else None
        if path is None or not path.is_file():
            fallback = run_directory / "Phaser" / f"mr_solution{suffix}"
            path = fallback if fallback.is_file() else None
        if path is None:
            raise TFZOverrideError(
                f"Phaser wrote no solution {name}; --tfz-override only accepts a low TFZ"
            )
    if tfz >= TFZ_THRESHOLD:
        return None
    return {
        "accepted_by": "--tfz-override",
        "mr_status": report.get("status"),
        "mr_message": report.get("message"),
        "best_tfz": float(tfz),
        "translation_tfz": phaser.get("translation_tfz"),
        "refined_tfz": phaser.get("refined_tfz"),
        "tfz_basis": phaser.get("tfz_basis"),   # None: recorded before both were kept
        "tfz_floor": TFZ_FLOOR,
        "tfz_threshold": TFZ_THRESHOLD,
        "accepted_utc": datetime.now(timezone.utc).isoformat(),
    }


def tfz_override_warning(acceptance: Mapping[str, object]) -> str:
    """One line for the terminal about an accepted low TFZ."""
    tfz = float(acceptance["best_tfz"])
    below = (f"below the {TFZ_FLOOR:.1f} floor" if tfz < TFZ_FLOOR
             else f"below the {TFZ_THRESHOLD:.1f} threshold")
    return (f"Warning: continuing from {acceptance.get('mr_status')} (best TFZ {tfz:.2f}, "
            f"{below}) because of --tfz-override")


_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)"
_TFZ_RE = re.compile(rf"\bTFZ\s*={{0,2}}\s*({_NUMBER})", re.I)
#: Phaser's solution annotation: "TFZ=" is the translation-function Z-score,
#: "TFZ==" the TFZ-equivalent after refinement, e.g.
#: "SOLU SET RFZ=3.9 TFZ=6.5 PAK=0 LLG=95 TFZ==7.8 LLG=96 TFZ==7.8".
_TRANSLATION_TFZ_RE = re.compile(rf"\bTFZ\s*=(?!=)\s*({_NUMBER})", re.I)
_REFINED_TFZ_RE = re.compile(rf"\bTFZ\s*==\s*({_NUMBER})", re.I)
_LLG_RE = re.compile(rf"\bLLG\s*={{0,2}}\s*({_NUMBER})", re.I)


def _phil_string(value: str | Path) -> str:
    """Return a double-quoted PHIL string, escaped the way libtbx.phil reads it.

    Inside quotes PHIL only decodes ``\\\\`` and ``\\"`` (libtbx
    ``quote_python_str``); any other backslash sequence is kept literally. JSON
    escaping therefore broke non-ASCII paths: ``\\u5b66`` reached Phaser as
    six literal characters and the file could not be found. Characters are
    written as they are; the parameter file is UTF-8.
    """
    text = str(value)
    if any(character in text for character in "\r\n\x00"):
        raise PhaserExecutionError(f"Path cannot be written to a PHIL file: {text!r}")
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _load_report(report_path: Path) -> dict[str, object]:
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PhaserExecutionError(f"Could not read frozen AutoMR report {report_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise PhaserExecutionError(f"Frozen AutoMR report is not a JSON object: {report_path}")
    if payload.get("workflow") != "automr" or payload.get("stage") != "preflight":
        raise PhaserExecutionError("Phaser requires a frozen AutoMR preflight report")
    return payload


def _integer_field(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise PhaserExecutionError(f"Frozen report has no valid {label}")
    return value


# Values written by the first (schema 1) version of the record.
_LEGACY_REQUESTED = {"auto": None, "intensities": INTENSITIES}


def frozen_phaser_data(report: Mapping[str, object]) -> dict[str, object]:
    """Return the preflight's data choice.

    Reports frozen before the record existed (and schema-1 ``auto`` records)
    have no labels: Phaser then keeps its own choice, as those runs did.
    """
    record = report.get("phaser_data")
    if record is None:
        return {"requested": None, "labin": None}
    if not isinstance(record, Mapping):
        raise PhaserExecutionError("Frozen phaser_data record is malformed")
    requested = record.get("requested")
    labin = record.get("labin")
    if record.get("schema_version") == 1 and requested in _LEGACY_REQUESTED:
        requested = _LEGACY_REQUESTED[requested]
        if requested is None and labin is None:
            return {**record, "requested": None}
    if (
        requested in {AMPLITUDES, INTENSITIES}
        and isinstance(labin, str)
        and re.fullmatch(r"[^\s,\"\\]+(?:,[^\s,\"\\]+)+", labin)
    ):
        return {**record, "requested": requested}
    raise PhaserExecutionError("Frozen phaser_data record is malformed")


def build_phaser_eff(report: Mapping[str, object], phaser_directory: Path) -> str:
    """Build one explicit ensemble, one NA composition, and one search request."""
    inputs = report.get("inputs")
    assessment = report.get("model_assessment")
    symmetry = report.get("symmetry")
    if not isinstance(inputs, Mapping) or not isinstance(assessment, Mapping):
        raise PhaserExecutionError("Frozen report is missing its input or model assessment")

    run = phaser_directory.parent
    reflections = resolve_artifact_path(inputs.get("reflections"), run)
    copied_model = resolve_artifact_path(assessment.get("copied_model"), run)
    if reflections is None:
        raise PhaserExecutionError("Frozen report reflections file is missing")
    if copied_model is None:
        raise PhaserExecutionError("Frozen, checksum-verified MR model is missing")
    model_sha256 = assessment.get("sha256")
    if not isinstance(model_sha256, str) or file_sha256(copied_model) != model_sha256:
        raise PhaserExecutionError("Frozen MR model does not match its recorded checksum")
    reflections_sha256 = inputs.get("reflections_sha256")
    if reflections_sha256 is not None and not isinstance(reflections_sha256, str):
        raise PhaserExecutionError("Frozen reflections checksum is malformed")
    if (
        isinstance(reflections_sha256, str)
        and file_sha256(reflections) != reflections_sha256
    ):
        raise PhaserExecutionError(
            "Frozen reflections do not match their recorded checksum"
        )

    residue_count = _integer_field(
        assessment.get("polymer_residue_count"), "polymer residue count"
    )
    copies = 1
    if isinstance(symmetry, Mapping) and "mr_copies" in symmetry:
        copies = _integer_field(symmetry.get("mr_copies"), "MR copy count")
    labin = frozen_phaser_data(report).get("labin")
    # Only reports frozen before phaser_data existed have no labin.
    data_lines = [f"  labin = {_phil_string(labin)}"] if labin else []

    # Do not set ``phaser.model`` as well as ``phaser.ensemble``.  The shortcut
    # can create an additional anonymous ensemble in some Phenix versions.
    return "\n".join([
        "phaser {",
        "  mode = MR_AUTO",
        f"  hklin = {_phil_string(str(reflections))}",
        *data_lines,
        "  chain_type = dna",
        "  composition {",
        "    chain {",
        "      chain_type = na",
        "      comp_type = nres",
        f"      nres = {residue_count}",
        f"      num = {copies}",
        "    }",
        "  }",
        "  ensemble {",
        "    model_id = nasolve_model",
        "    use_hetatm = True",
        "    coordinates {",
        f"      pdb = {_phil_string(str(copied_model))}",
        "      rmsd = 1.0",
        "    }",
        "  }",
        "  search {",
        "    ensembles = nasolve_model",
        f"    copies = {copies}",
        "  }",
        f"  output_dir = {_phil_string(phaser_directory.resolve())}",
        "  keywords {",
        "    general {",
        '      root = "PHASER"',
        "    }",
        "  }",
        "}",
        "",
    ])


@dataclass(frozen=True)
class TFZScores:
    """The best Phaser solution's TFZ values.

    ``tfz`` decides pass/fail: the TFZ after refinement (``TFZ==``, the last on
    the solution line) when Phaser reports one, otherwise the translation-search
    TFZ (``TFZ=``, the last on the line: the last placed copy). ``llg`` is the
    last LLG on the same line.
    """

    tfz: float | None
    llg: float | None
    translation_tfz: float | None
    refined_tfz: float | None
    basis: str | None   # "refined", "translation" or None


def _last(pattern: re.Pattern[str], line: str) -> float | None:
    values = pattern.findall(line)
    return float(values[-1]) if values else None


def parse_tfz_scores(log_text: str) -> TFZScores:
    """Return the TFZ values of the best final solution (highest deciding TFZ)."""
    lines = log_text.splitlines()
    candidates = [line for line in lines
                  if re.search(r"\bSOLU\s+SET\b", line, re.I) and _TFZ_RE.search(line)]
    candidates = candidates or [line for line in lines if _TFZ_RE.search(line)]
    best: TFZScores | None = None
    for line in candidates:
        refined = _last(_REFINED_TFZ_RE, line)
        translation = _last(_TRANSLATION_TFZ_RE, line)
        deciding = refined if refined is not None else translation
        if deciding is None:
            continue
        scores = TFZScores(
            tfz=deciding,
            llg=_last(_LLG_RE, line),
            translation_tfz=translation,
            refined_tfz=refined,
            basis="refined" if refined is not None else "translation",
        )
        if best is None or scores.tfz > best.tfz:
            best = scores
    return best or TFZScores(None, None, None, None, None)


def parse_best_tfz(log_text: str) -> tuple[float | None, float | None]:
    """Return the deciding TFZ of the best final solution and its LLG."""
    scores = parse_tfz_scores(log_text)
    return scores.tfz, scores.llg


def describe_tfz(scores: TFZScores) -> str:
    """One line with both TFZ values, saying which decides pass/fail."""
    def number(value: float | None) -> str:
        return f"{value:.2f}" if value is not None else "not reported"

    if scores.basis == "refined":
        return (f"TFZ {number(scores.refined_tfz)} after refinement (decides pass/fail); "
                f"{number(scores.translation_tfz)} from the translation search")
    if scores.basis == "translation":
        return (f"TFZ {number(scores.translation_tfz)} from the translation search "
                "(decides pass/fail; Phaser reported no TFZ after refinement)")
    return "TFZ not reported"


_BASIS_WORDS = {"refined": "after refinement", "translation": "translation search"}


def _basis_note(scores: TFZScores) -> str:
    if scores.basis == "refined" and scores.translation_tfz is not None:
        return f" (after refinement; translation search {scores.translation_tfz:.2f})"
    if scores.basis == "refined":
        return " (after refinement)"
    if scores.basis == "translation":
        return " (translation search; no TFZ after refinement reported)"
    return ""


def _solution_sort_key(path: Path) -> tuple[int, str]:
    match = re.search(r"\.(\d+)\.(?:pdb|mtz)$", path.name, re.I)
    return (int(match.group(1)) if match else 10**9, path.name.casefold())


def _find_solution(phaser_directory: Path, suffix: str) -> Path | None:
    candidates = sorted(
        (
            path for path in phaser_directory.glob(f"PHASER*{suffix}")
            if path.is_file() and not path.name.startswith("mr_solution")
        ),
        key=_solution_sort_key,
    )
    return candidates[0] if candidates else None


def _write_report(path: Path, payload: Mapping[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _append_run_log(path: Path, lines: list[str]) -> None:
    try:
        existing = path.read_text(encoding="utf-8")
    except OSError:
        existing = ""
    existing = existing.replace("Phaser executed: no", "Phaser executed: yes", 1)
    existing = existing.replace(
        "Next layer: use this frozen run input to prepare and execute Phenix Phaser.",
        "Phaser execution completed; details follow.",
        1,
    )
    if existing and not existing.endswith("\n"):
        existing += "\n"
    path.write_text(existing + "\n".join(lines) + "\n", encoding="utf-8")


def execute_phaser(
    report_path: Path,
    phaser_executable: Path,
    environment: Mapping[str, str] | None = None,
    phenix_version: str = "unknown",
) -> PhaserExecutionResult:
    """Run Phaser once, preserve its raw output, and apply the AutoMR TFZ gate."""
    frozen_report_path = report_path.expanduser().resolve()
    report = _load_report(frozen_report_path)
    run_directory = frozen_report_path.parent
    phaser_directory = run_directory / "Phaser"
    try:
        phaser_directory.mkdir()
    except FileExistsError as exc:
        raise PhaserExecutionError(
            f"Phaser directory already exists; refusing to overwrite {phaser_directory}"
        ) from exc

    parameter_path = phaser_directory / "phaser.eff"
    log_path = phaser_directory / "phaser.log"
    parameter_path.write_text(build_phaser_eff(report, phaser_directory), encoding="utf-8")
    executable = phaser_executable.expanduser().resolve()
    command = [str(executable), str(parameter_path)]
    try:
        completed = subprocess.run(
            command,
            cwd=phaser_directory,
            env=dict(environment) if environment is not None else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        output = completed.stdout
        return_code = completed.returncode
    except OSError as exc:
        output = f"NASolve could not launch Phaser: {exc}\n"
        return_code = 127
    log_path.write_text(output, encoding="utf-8")

    scores = parse_tfz_scores(output)
    tfz, llg = scores.tfz, scores.llg
    phaser_data = frozen_phaser_data(report)
    labin_echo = reported_labin(output)
    used_kind = reported_data_kind(labin_echo)
    requested_kind = phaser_data.get("requested")
    data_mismatch = (
        requested_kind is not None
        and used_kind is not None
        and used_kind != requested_kind
    )
    raw_pdb = _find_solution(phaser_directory, ".pdb")
    raw_mtz = _find_solution(phaser_directory, ".mtz")
    solution_pdb: Path | None = None
    solution_mtz: Path | None = None
    if raw_pdb:
        solution_pdb = phaser_directory / "mr_solution.pdb"
        shutil.copyfile(raw_pdb, solution_pdb)
    if raw_mtz:
        solution_mtz = phaser_directory / "mr_solution.mtz"
        shutil.copyfile(raw_mtz, solution_mtz)

    outputs_complete = solution_pdb is not None and solution_mtz is not None
    if return_code != 0:
        status = "MR_FAILED"
        message = f"Phaser exited with status {return_code}; inspect phaser.log"
    elif data_mismatch:
        status = "MR_FAILED"
        message = (
            f"{describe(str(requested_kind)).capitalize()} ({phaser_data.get('labin')}) "
            f"were requested but Phaser reported LABIN {labin_echo}; inspect phaser.log"
        )
    elif tfz is None:
        status = "MR_FAILED"
        message = "Phaser did not report a TFZ; inspect phaser.log"
    elif tfz < TFZ_FLOOR:
        status = "MR_FAILED"
        message = (f"Best TFZ {tfz:.2f}{_basis_note(scores)} is below the "
                   f"{TFZ_FLOOR:.1f} acceptance floor")
    elif not outputs_complete:
        status = "MR_FAILED"
        message = "Phaser reported a candidate solution but did not produce both PDB and MTZ outputs"
    elif tfz < TFZ_THRESHOLD:
        status = "MR_REVIEW"
        message = f"Best TFZ {tfz:.2f}{_basis_note(scores)} is a red flag and requires user review"
    else:
        status = "MR_SUCCESS"
        message = (f"Best TFZ {tfz:.2f}{_basis_note(scores)} passes the "
                   f"{TFZ_THRESHOLD:.1f} AutoMR threshold")

    previous_status = report.get("status")
    previous_message = report.get("message")
    report["stage"] = "phaser"
    report["status"] = status
    report["message"] = message
    report["updated_utc"] = datetime.now(timezone.utc).isoformat()
    report["preflight"] = {"status": previous_status, "message": previous_message}
    report["execution"] = {
        "phaser_ran": True,
        "phaser": {
            "executable": str(executable),
            "phenix_version": phenix_version,
            "command": command,
            "return_code": return_code,
            "parameters": str(parameter_path),
            "log": str(log_path),
            "ensemble_count": 1,
            "ensemble_id": "nasolve_model",
            "preserve_heteroatoms": True,
            "model_rmsd": 1.0,
            "composition_nres": report["model_assessment"]["polymer_residue_count"],
            "reflection_data": {
                "requested": phaser_data.get("requested"),
                "labin": phaser_data.get("labin"),
                "phaser_reported_labin": labin_echo,
                "phaser_reported_kind": used_kind,
            },
            "search_copies": (
                report["symmetry"].get("mr_copies", 1)
                if isinstance(report.get("symmetry"), Mapping) else 1
            ),
            "best_tfz": tfz,   # the deciding TFZ (tfz_basis)
            "best_llg": llg,
            "translation_tfz": scores.translation_tfz,
            "refined_tfz": scores.refined_tfz,
            "tfz_basis": scores.basis,
            "raw_solution_pdb": str(raw_pdb) if raw_pdb else None,
            "raw_solution_mtz": str(raw_mtz) if raw_mtz else None,
            "solution_pdb": str(solution_pdb) if solution_pdb else None,
            "solution_mtz": str(solution_mtz) if solution_mtz else None,
        },
    }
    _write_report(frozen_report_path, report)
    _append_run_log(
        run_directory / "automr.log",
        [
            "Phaser execution",
            f"Status: {status}",
            f"Phenix version: {phenix_version}",
            "Ensembles: 1 (nasolve_model)",
            "Preserve heteroatoms: True",
            "Model RMSD: 1.0 A",
            "Reflection data: "
            + (
                f"{requested_kind} - {describe(str(requested_kind))} "
                f"({phaser_data.get('labin')})"
                if phaser_data.get("labin") else "Phenix default (legacy report)"
            )
            + (f"; Phaser reported LABIN {labin_echo}" if labin_echo else ""),
            f"Composition residues: {report['model_assessment']['polymer_residue_count']}",
            f"Search copies: {report['execution']['phaser']['search_copies']}",
            f"Best TFZ: {tfz if tfz is not None else 'not reported'}"
            + (f" ({_BASIS_WORDS[scores.basis]}; decides pass/fail)" if scores.basis else ""),
            "TFZ after refinement (TFZ==): "
            + (f"{scores.refined_tfz}" if scores.refined_tfz is not None else "not reported"),
            "TFZ from the translation search (TFZ=): "
            + (f"{scores.translation_tfz}" if scores.translation_tfz is not None else "not reported"),
            f"Best LLG: {llg if llg is not None else 'not reported'}",
            f"Result: {message}",
            f"Phaser log: {log_path}",
            "",
        ],
    )
    return PhaserExecutionResult(
        status=status,
        message=message,
        run_directory=run_directory,
        phaser_directory=phaser_directory,
        report_path=frozen_report_path,
        log_path=log_path,
        parameter_path=parameter_path,
        tfz=tfz,
        llg=llg,
        solution_pdb=solution_pdb,
        solution_mtz=solution_mtz,
        translation_tfz=scores.translation_tfz,
        refined_tfz=scores.refined_tfz,
        tfz_basis=scores.basis,
    )


__all__ = [
    "TFZOverrideError",
    "PhaserExecutionError",
    "PhaserExecutionResult",
    "TFZScores",
    "TFZ_FLOOR",
    "TFZ_THRESHOLD",
    "build_phaser_eff",
    "execute_phaser",
    "frozen_phaser_data",
    "check_tfz_override",
    "tfz_override_warning",
    "describe_tfz",
    "parse_best_tfz",
    "parse_tfz_scores",
]
