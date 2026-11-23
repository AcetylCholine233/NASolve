"""Numbered PostMR results of one AutoMR run, and switching between them.

An AutoMR run holds one *current* PostMR result: ``PostMR/`` and everything
made from it (``AutoSol/``, ``AutoRefine/`` with its rounds and checkpoints,
``RefineDoctor/``, ``CootGUI/{postmr,autosol,autorefine}``), with the matching
run-report records (``postmr``, ``autosol``, ``autorefine``,
``refine_doctor``) and the run's stage, status and message.

``nasolve postmr RUN --force-overwrite`` keeps the current result before
running PostMR again (user's request, 2026-11-23; it first deleted it): its
folders are *moved*, unchanged, to ``PostMRHistory/postmr-NNN/`` (same layout)
and its report records saved in ``PostMRHistory/postmr-NNN/result.json``.
``nasolve postmr RUN --switch postmr-NNN`` puts the current result away the
same way and moves the chosen one back. ``--history`` lists them.

Results are moved, never copied or edited, and always used from the run's own
folders: every path a result recorded (absolute paths, run-anchored references
with checksums, the paths inside phenix.refine's .eff files) is valid again
once it is current. While a result is in the history its recorded paths point
at the current result's folders, so open it by switching to it. The Phaser
solution, metal searches and merges belong to the run and are never moved.
"""

from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

HISTORY_FOLDER = "PostMRHistory"
#: What a PostMR result consists of, relative to the AutoMR run.
RESULT_FOLDERS = ("PostMR", "AutoSol", "AutoRefine", "RefineDoctor",
                  "CootGUI/postmr", "CootGUI/autosol", "CootGUI/autorefine")
RESULT_RECORDS = ("postmr", "autosol", "autorefine", "refine_doctor")
RESULT_FIELDS = ("stage", "status", "message", "updated_utc")
RESULT_FILE = "result.json"
_ID = re.compile(r"postmr-(\d{3,})")


class PostMRHistoryError(RuntimeError):
    """A PostMR result cannot be kept, listed or switched to."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def current_id(report: Mapping[str, object]) -> str | None:
    """The number of the current result (``report["postmr_result"]["id"]``)."""
    block = report.get("postmr_result")
    value = block.get("id") if isinstance(block, Mapping) else None
    return value if isinstance(value, str) and _ID.fullmatch(value) else None


def kept_ids(run: Path) -> list[str]:
    """Results kept in ``PostMRHistory/``, in number order."""
    folder = run / HISTORY_FOLDER
    if not folder.is_dir():
        return []
    found = [path.name for path in folder.iterdir()
             if _ID.fullmatch(path.name) and (path / RESULT_FILE).is_file()]
    return sorted(found, key=lambda name: int(_ID.fullmatch(name).group(1)))


def next_id(run: Path, report: Mapping[str, object]) -> str:
    """A number no result of this run has used."""
    folder = run / HISTORY_FOLDER
    used = [int(match.group(1)) for path in (folder.iterdir() if folder.is_dir() else ())
            if (match := _ID.fullmatch(path.name))]
    if (active := current_id(report)) is not None:
        used.append(int(_ID.fullmatch(active).group(1)))
    return f"postmr-{max(used, default=0) + 1:03d}"


def has_current(run: Path, report: Mapping[str, object]) -> bool:
    return any((run / folder).exists() or (run / folder).is_symlink() for folder in RESULT_FOLDERS) \
        or any(key in report for key in RESULT_RECORDS)


def _move(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    os.rename(source, destination)


def _summary(root: Path, records: Mapping[str, object]) -> dict[str, object]:
    """What ``--history`` shows for a result whose folders are under ``root``."""
    postmr = records.get("postmr") if isinstance(records.get("postmr"), Mapping) else {}
    merge = postmr.get("metal_merge") if isinstance(postmr.get("metal_merge"), Mapping) else {}
    rounds = sorted(path.name for path in (root / "AutoRefine").glob("round_*") if path.is_dir()) \
        if (root / "AutoRefine").is_dir() else []
    current_checkpoint = metrics = None
    try:
        registry = json.loads((root / "AutoRefine" / "checkpoints.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        registry = None
    if isinstance(registry, Mapping):
        current_checkpoint = registry.get("current")
        for item in registry.get("checkpoints") or []:
            if isinstance(item, Mapping) and item.get("id") == current_checkpoint:
                values = item.get("metrics")
                if isinstance(values, Mapping):
                    metrics = {key: values.get(key) for key in ("r_work", "r_free")}
    model = postmr.get("input_model")
    return {
        "postmr_created_utc": postmr.get("created_utc"),
        "postmr_status": postmr.get("status"),
        "starting_model": Path(model).name if isinstance(model, str) else None,
        "merged_model": merge.get("name"),
        "autorefine_rounds": len(rounds),
        "current_checkpoint": current_checkpoint,
        "current_metrics": metrics,
    }


def keep_current(run: Path, report: dict[str, object], *, reason: str) -> dict[str, object] | None:
    """Move the current result to ``PostMRHistory/postmr-NNN/`` and take its
    records out of ``report`` (not written here). Returns the kept record, or
    None when the run has no current result."""
    if not has_current(run, report):
        return None
    result_id = current_id(report) or next_id(run, report)
    target = run / HISTORY_FOLDER / result_id
    if target.exists():
        raise PostMRHistoryError(f"{target} already exists; refusing to overwrite it")
    present = [folder for folder in RESULT_FOLDERS
               if (run / folder).exists() or (run / folder).is_symlink()]
    records = {key: report[key] for key in RESULT_RECORDS if key in report}
    entry = {
        "schema_version": 1,
        "id": result_id,
        "kept_utc": _now(),
        "reason": reason,
        "folders": present,
        "records": records,
        "fields": {key: report.get(key) for key in RESULT_FIELDS},
        "postmr_result": report.get("postmr_result"),
        "summary": _summary(run, records),
    }
    moved: list[str] = []
    try:
        target.mkdir(parents=True)
        for folder in present:
            _move(run / folder, target / folder)
            moved.append(folder)
        _write_json(target / RESULT_FILE, entry)
    except OSError as exc:
        for folder in reversed(moved):   # put everything back
            _move(target / folder, run / folder)
        if target.is_dir():
            _remove_empty_tree(target)
        raise PostMRHistoryError(f"Could not keep the current PostMR result as {result_id}: {exc}") from exc
    for key in (*RESULT_RECORDS, "postmr_result"):
        report.pop(key, None)
    return entry


def _remove_empty_tree(folder: Path) -> None:
    for path in sorted(folder.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_dir() and not any(path.iterdir()):
            path.rmdir()
    if folder.is_dir() and not any(folder.iterdir()):
        folder.rmdir()


def _restore(run: Path, report: dict[str, object], result_id: str) -> dict[str, object]:
    target = run / HISTORY_FOLDER / result_id
    try:
        entry = json.loads((target / RESULT_FILE).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise PostMRHistoryError(f"Cannot read {target / RESULT_FILE}: {exc}") from exc
    folders = [folder for folder in entry.get("folders") or [] if folder in RESULT_FOLDERS]
    in_the_way = [folder for folder in folders if (run / folder).exists()]
    if in_the_way:
        raise PostMRHistoryError(f"Cannot restore {result_id}: {', '.join(in_the_way)} already exist")
    moved: list[str] = []
    try:
        for folder in folders:
            _move(target / folder, run / folder)
            moved.append(folder)
    except OSError as exc:
        for folder in reversed(moved):
            _move(run / folder, target / folder)
        raise PostMRHistoryError(f"Could not restore {result_id}: {exc}") from exc
    records = entry.get("records") if isinstance(entry.get("records"), Mapping) else {}
    for key in RESULT_RECORDS:
        report.pop(key, None)
    report.update({key: value for key, value in records.items() if key in RESULT_RECORDS})
    fields = entry.get("fields") if isinstance(entry.get("fields"), Mapping) else {}
    report.update({key: fields.get(key) for key in RESULT_FIELDS if key in fields})
    previous = entry.get("postmr_result") if isinstance(entry.get("postmr_result"), Mapping) else {}
    report["postmr_result"] = {**previous, "id": result_id, "restored_utc": _now()}
    (target / RESULT_FILE).unlink()
    _remove_empty_tree(target)
    return entry


def switch_to(run_directory: Path, result_id: str) -> tuple[dict[str, object] | None, dict[str, object]]:
    """Make the kept result ``result_id`` current again; the current result
    (if any) is kept first. Returns (kept current, restored) records."""
    run = run_directory.expanduser().resolve()
    report_path = run / "report.json"
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise PostMRHistoryError(f"Could not read AutoMR report {report_path}: {exc}") from exc
    if not isinstance(report, dict) or report.get("workflow") != "automr":
        raise PostMRHistoryError(f"{report_path} is not an AutoMR run report")
    if result_id == current_id(report) and has_current(run, report):
        raise PostMRHistoryError(f"{result_id} is already the current PostMR result")
    available = kept_ids(run)
    if result_id not in available:
        raise PostMRHistoryError(
            f"No kept PostMR result {result_id} in {run / HISTORY_FOLDER}; kept: "
            + (", ".join(available) or "none")
        )
    kept = keep_current(run, report, reason=f"nasolve postmr --switch {result_id}")
    try:
        restored = _restore(run, report, result_id)
    except PostMRHistoryError:
        if kept is not None:   # put the previous current result back
            _restore(run, report, str(kept["id"]))
            _write_json(report_path, report)
        raise
    switches = list(report.get("postmr_switches") or [])
    switches.append({"switched_utc": _now(), "to": result_id,
                     "kept": kept["id"] if kept is not None else None})
    report["postmr_switches"] = switches
    _write_json(report_path, report)
    return kept, restored


def history(run_directory: Path) -> list[dict[str, object]]:
    """The run's PostMR results, the current one marked ``current``."""
    run = run_directory.expanduser().resolve()
    try:
        report = json.loads((run / "report.json").read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise PostMRHistoryError(f"Could not read AutoMR report {run / 'report.json'}: {exc}") from exc
    rows: list[dict[str, object]] = []
    for result_id in kept_ids(run):
        entry = json.loads((run / HISTORY_FOLDER / result_id / RESULT_FILE).read_text(encoding="utf-8"))
        records = entry.get("records") if isinstance(entry.get("records"), Mapping) else {}
        rows.append({"id": result_id, "current": False, "kept_utc": entry.get("kept_utc"),
                     "fields": entry.get("fields"),
                     **_summary(run / HISTORY_FOLDER / result_id, records)})
    if isinstance(report, Mapping) and has_current(run, report):
        rows.append({"id": current_id(report), "current": True, "kept_utc": None,
                     "fields": {key: report.get(key) for key in RESULT_FIELDS},
                     **_summary(run, report)})
    return sorted(rows, key=lambda row: (row["id"] is None, row["id"] or ""))


__all__ = [
    "HISTORY_FOLDER", "PostMRHistoryError", "RESULT_FOLDERS", "RESULT_RECORDS",
    "current_id", "has_current", "history", "keep_current", "kept_ids", "next_id", "switch_to",
]
