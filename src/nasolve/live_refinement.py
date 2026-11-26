"""Live snapshots of a running phenix.refine, for looking at it in Coot.

phenix.refine run from the command line writes its model only at the end.
Phenix's own GUI shows a refinement as it goes through a hook in phenix.refine
itself (Phenix 2.2.1, ``phenix/refinement/driver.py`` and ``runtime.py``):
when the environment variable ``PHENIX_GUI_ENVIRONMENT`` is set, it saves a
snapshot, ``.status.pkl`` in its working directory, at the start, after every
macro-cycle and at the end. The snapshot holds the current model (coordinates,
B, occupancy) and 2mFo-DFc and mFo-DFc map coefficients (filled). The variable
changes nothing else in phenix.refine except the wording of one R-free-flag
error message (``iotbx/extract_xtal_data.py``).

AutoRefine sets the variable and, while phenix.refine runs, turns each snapshot
into ``round_NNN/live/snapshot_NNN.pdb`` and ``snapshot_NNN_maps.mtz``
(2FOFCWT/PH2FOFCWT, FOFCWT/PHFOFCWT) with phenix.python, and lists them in
``round_NNN/live/status.json`` for ``nasolve show RUN --live``. The snapshot
is read from a copy; conversion runs beside the refinement and its failures
never change the refinement or its result.
"""

from __future__ import annotations

import json
import math
import os
import queue
import re
import shutil
import subprocess
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Mapping, Sequence

#: The hook phenix.refine checks (Phenix 2.2.1 ``phenix/refinement/driver.py``).
GUI_ENVIRONMENT_VARIABLE = "PHENIX_GUI_ENVIRONMENT"
#: What phenix.refine writes in its working directory when the variable is set.
PHENIX_STATUS_FILE = ".status.pkl"
LIVE_DIRECTORY = "live"
STATUS_FILE = "status.json"
SCHEMA_VERSION = 1
MAP_LABELS = (("2FOFCWT", "PH2FOFCWT"), ("FOFCWT", "PHFOFCWT"))
METHOD = (
    "PHENIX_GUI_ENVIRONMENT=1: phenix.refine saves .status.pkl at the start, after "
    "each macro-cycle and at the end (the hook Phenix's GUI uses)"
)
#: Seconds between looks at .status.pkl; a snapshot is copied once it has
#: stayed the same for one look.
POLL_SECONDS = 0.5
CONVERT_TIMEOUT = 120
#: After phenix.refine ends, how long to wait for snapshots still converting.
FINISH_TIMEOUT = 180
#: Conversion stops after this many failures in a row.
MAX_FAILURES = 3

_MARKER = "NASOLVE_LIVE_SNAPSHOT:"
#: Run by phenix.python: snapshot copy -> model PDB + map coefficients MTZ.
CONVERTER = r'''
import json, math, os, sys
from libtbx import easy_pickle
import iotbx.map_tools
source, model_path, maps_path = sys.argv[1:4]
status = easy_pickle.load(source)
symmetry = status.xray_structure.crystal_symmetry()
status.pdb_hierarchy.write_pdb_file(file_name=model_path + ".part", crystal_symmetry=symmetry)
os.replace(model_path + ".part", model_path)
iotbx.map_tools.write_map_coeffs(status.two_fofc_coeffs, status.fofc_coeffs,
                                 file_name=maps_path + ".part")
os.replace(maps_path + ".part", maps_path)
def number(name):
    try:
        value = float(getattr(status, name, None))
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None
print("NASOLVE_LIVE_SNAPSHOT:" + json.dumps({
    "r_work": number("r_work"), "r_free": number("r_free"),
    "bonds": number("b_mean"), "angles": number("a_mean")}))
'''


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _relative(path: Path, base: Path) -> str:
    """``path`` relative to ``base`` when it can be written that way, else absolute."""
    try:
        return os.path.relpath(path.resolve(), base.resolve())
    except ValueError:   # another drive (Windows)
        return str(path.resolve())


def snapshot_label(number: int, macro_cycles: int) -> tuple[str, int]:
    """(label, macro-cycles completed) of the ``number``-th snapshot.

    phenix.refine saves one before macro-cycle 1, one after each macro-cycle
    and one at the end.
    """
    if number <= 1:
        return "start", 0
    if number <= macro_cycles + 1:
        return f"macro-cycle {number - 1} of {macro_cycles}", number - 1
    return "final", macro_cycles


def format_snapshot(checkpoint_id: str, snapshot: Mapping[str, object]) -> str:
    def value(key: str) -> str:
        number = snapshot.get(key)
        return f"{number:.4f}" if isinstance(number, float) else "?"
    return (f"  {checkpoint_id} {snapshot.get('label')}: "
            f"R-work {value('r_work')}, R-free {value('r_free')}")


class LiveSnapshots:
    """Watches one phenix.refine run and converts its snapshots.

    ``start`` before phenix.refine is launched, ``poll`` while it runs (from the
    thread waiting on it), ``finish`` when it has ended. Every method swallows
    its own errors into ``warnings``: the live view must never break a round.
    """

    def __init__(
        self,
        round_directory: Path,
        phenix_python: Path | None,
        environment: Mapping[str, str] | None,
        *,
        checkpoint_id: str,
        parent_checkpoint: str | None,
        start_model: Path,
        dictionaries: Sequence[Path] = (),
        macro_cycles: int,
        report: Callable[[str], None] | None = None,
    ) -> None:
        self.round_directory = Path(round_directory)
        self.directory = self.round_directory / LIVE_DIRECTORY
        self.status_path = self.directory / STATUS_FILE
        self.source = self.round_directory / PHENIX_STATUS_FILE
        self.phenix_python = phenix_python
        self.environment = dict(environment) if environment is not None else None
        self.checkpoint_id = checkpoint_id
        self.parent_checkpoint = parent_checkpoint
        self.start_model = Path(start_model)
        self.dictionaries = [Path(path) for path in dictionaries]
        self.macro_cycles = macro_cycles
        self.report = report
        self.reason: str | None = None
        if phenix_python is None or not Path(phenix_python).is_file():
            self.reason = f"phenix.python not found next to phenix.refine ({phenix_python})"
        self.enabled = self.reason is None
        self.state = "not-started"
        self.snapshots: list[dict[str, object]] = []
        self.warnings: list[str] = []
        self.copies = 0
        self.failures = 0
        self.converting = self.enabled
        self._seen: tuple[int, int] | None = None
        self._pending: tuple[int, int] | None = None
        self._lock = threading.Lock()
        self._queue: queue.Queue = queue.Queue()
        self._worker: threading.Thread | None = None
        self.started_utc: str | None = None

    # -- phenix.refine's environment ------------------------------------
    def refine_environment(self, base: Mapping[str, str] | None) -> dict[str, str] | None:
        """phenix.refine's environment: ``base`` plus the snapshot hook when enabled."""
        if not self.enabled:
            return dict(base) if base is not None else None
        environment = dict(os.environ if base is None else base)
        environment[GUI_ENVIRONMENT_VARIABLE] = "1"
        return environment

    # -- life cycle -----------------------------------------------------
    def start(self) -> None:
        self.started_utc = _now()
        if not self.enabled:
            self.state = "unavailable"
            return
        try:
            self.directory.mkdir()
            self._seen = self._stat()   # an old snapshot (none expected) is not new
            self.state = "running"
            self._write_status()
            self._worker = threading.Thread(target=self._work, name="nasolve-live", daemon=True)
            self._worker.start()
        except Exception as exc:   # pragma: no cover - defensive
            self._disable(f"live view not started: {exc}")

    def poll(self) -> None:
        if not self.enabled or self.state != "running":
            return
        try:
            self._poll(final=False)
        except Exception as exc:   # pragma: no cover - defensive
            self._warn(f"could not read the snapshot: {exc}")

    def finish(self, returncode: int | None, *, interrupted: bool = False) -> dict[str, object]:
        """Collect the last snapshot, wait for conversions, record the end."""
        if self.enabled and self.state == "running":
            try:
                if not interrupted:
                    self._poll(final=True)
                self._queue.put(None)
                if self._worker is not None and not interrupted:
                    self._worker.join(timeout=FINISH_TIMEOUT)
                    if self._worker.is_alive():
                        self._warn("stopped waiting for snapshots still being converted")
            except Exception as exc:   # pragma: no cover - defensive
                self._warn(f"could not finish the live view: {exc}")
            self.state = ("interrupted" if interrupted
                          else "finished" if returncode == 0 else "failed")
            self._write_status()
        return self.record()

    def record(self) -> dict[str, object]:
        """The round report's ``live`` block."""
        with self._lock:
            return {
                "enabled": self.enabled,
                "reason": self.reason,
                "method": METHOD if self.enabled else None,
                "status_file": (_relative(self.status_path, self.round_directory)
                                if self.enabled else None),
                "state": self.state,
                "snapshots": len(self.snapshots),
                "warnings": list(self.warnings),
            }

    # -- internals ------------------------------------------------------
    def _stat(self) -> tuple[int, int] | None:
        try:
            stat = self.source.stat()
        except OSError:
            return None
        return stat.st_mtime_ns, stat.st_size

    def _poll(self, *, final: bool) -> None:
        current = self._stat()
        if current is None or current == self._seen:
            self._pending = None
            return
        if current != self._pending and not final:
            self._pending = current      # changed: copy once it has settled
            return
        copy = self.directory / f".snapshot_{self.copies + 1:03d}.pkl"
        shutil.copyfile(self.source, copy)
        if self._stat() != current:      # rewritten while copying
            copy.unlink()
            self._pending = None
            return
        self._seen = current
        self._pending = None
        self.copies += 1
        label, completed = snapshot_label(self.copies, self.macro_cycles)
        self._queue.put({"number": self.copies, "label": label, "macro_cycle": completed,
                         "copy": copy, "taken_utc": _now()})

    def _work(self) -> None:
        while True:
            item = self._queue.get()
            if item is None:
                return
            if not self.converting:
                continue
            try:
                snapshot = self._convert(item)
            except Exception as exc:
                self.failures += 1
                self._warn(f"snapshot {item['number']} ({item['label']}) not converted: {exc}")
                if self.failures >= MAX_FAILURES:
                    self.converting = False
                    self._warn(f"live view stopped after {MAX_FAILURES} failed conversions; "
                               "the refinement itself is not affected")
                self._write_status()
                continue
            self.failures = 0
            with self._lock:
                self.snapshots.append(snapshot)
            self._write_status()
            try:
                item["copy"].unlink()
            except OSError:
                pass
            if self.report is not None:
                try:
                    self.report(format_snapshot(self.checkpoint_id, snapshot))
                except Exception:   # pragma: no cover - printing must not stop the watcher
                    pass

    def _convert(self, item: Mapping[str, object]) -> dict[str, object]:
        stem = f"snapshot_{int(item['number']):03d}"
        model = self.directory / f"{stem}.pdb"
        maps = self.directory / f"{stem}_maps.mtz"
        try:
            completed = subprocess.run(
                [str(self.phenix_python), "-c", CONVERTER, str(item["copy"]), str(model), str(maps)],
                env=self.environment, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=False, timeout=CONVERT_TIMEOUT,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(f"phenix.python took longer than {CONVERT_TIMEOUT} s") from exc
        matches = re.findall(rf"^{re.escape(_MARKER)}(\{{.*\}})$", completed.stdout, re.M)
        if completed.returncode or len(matches) != 1 or not model.is_file() or not maps.is_file():
            tail = " ".join(completed.stdout.strip().splitlines()[-1:])[:200]
            raise RuntimeError(f"phenix.python exited {completed.returncode}" + (f": {tail}" if tail else ""))
        values = json.loads(matches[0])
        numbers = {
            key: (float(values[key]) if isinstance(values.get(key), (int, float))
                  and not isinstance(values.get(key), bool) and math.isfinite(values[key]) else None)
            for key in ("r_work", "r_free", "bonds", "angles")
        }
        return {
            "number": item["number"], "label": item["label"], "macro_cycle": item["macro_cycle"],
            "model": model.name, "maps": maps.name, **numbers,
            "taken_utc": item["taken_utc"], "written_utc": _now(),
        }

    def _warn(self, message: str) -> None:
        with self._lock:
            self.warnings.append(message)
        if self.report is not None:
            try:
                self.report(f"  {self.checkpoint_id} live view: {message}")
            except Exception:   # pragma: no cover
                pass

    def _disable(self, message: str) -> None:
        self.enabled = False
        self.reason = message
        self.state = "unavailable"

    def _write_status(self) -> None:
        with self._lock:
            payload = {
                "schema_version": SCHEMA_VERSION,
                "kind": "nasolve-autorefine-live",
                "checkpoint": self.checkpoint_id,
                "parent_checkpoint": self.parent_checkpoint,
                "round": self.round_directory.name,
                "state": self.state,
                "pid": os.getpid(),
                "macro_cycles": self.macro_cycles,
                # Relative to this live directory when possible.
                "start_model": _relative(self.start_model, self.directory),
                "dictionaries": [_relative(path, self.directory) for path in self.dictionaries],
                "map_labels": [list(pair) for pair in MAP_LABELS],
                "method": METHOD,
                "started_utc": self.started_utc,
                "updated_utc": _now(),
                "snapshots": list(self.snapshots),
                "warnings": list(self.warnings),
            }
            temporary = self.status_path.with_name(self.status_path.name + ".part")
            temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            os.replace(temporary, self.status_path)


def run_with_snapshots(
    command: Sequence[str],
    *,
    cwd: Path,
    environment: Mapping[str, str] | None,
    log_handle,
    watcher: LiveSnapshots,
) -> int:
    """Run phenix.refine like ``subprocess.run`` while ``watcher`` collects snapshots."""
    watcher.start()
    try:
        process = subprocess.Popen(
            list(command), cwd=cwd, env=watcher.refine_environment(environment), text=True,
            stdout=log_handle, stderr=subprocess.STDOUT,
        )
    except OSError:
        watcher.finish(127)
        raise
    try:
        while True:
            try:
                returncode = process.wait(timeout=POLL_SECONDS)
                break
            except subprocess.TimeoutExpired:
                watcher.poll()
    except BaseException:
        # Ctrl+C (or another error): stop phenix.refine as subprocess.run would.
        try:
            process.wait(timeout=0.25)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        watcher.finish(process.returncode, interrupted=True)
        raise
    watcher.finish(returncode)
    return returncode


# ---------------------------------------------------------------------------
# Reading a live view (nasolve show RUN --live)


class LiveViewError(RuntimeError):
    """No live snapshot can be shown."""


def _process_alive(pid: object) -> bool:
    if not isinstance(pid, int) or isinstance(pid, bool) or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def _plain_name(value: object, what: str) -> str:
    if not isinstance(value, str) or not value or Path(value).name != value or value in {".", ".."}:
        raise LiveViewError(f"Live status names an invalid {what}: {value!r}")
    return value


def read_live_status(status_path: Path) -> dict[str, object]:
    try:
        payload = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise LiveViewError(f"Cannot read {status_path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("kind") != "nasolve-autorefine-live" \
            or payload.get("schema_version") != SCHEMA_VERSION:
        raise LiveViewError(f"{status_path} is not a NASolve live status file")
    if payload.get("state") == "running" and not _process_alive(payload.get("pid")):
        payload["state"] = "stopped"   # NASolve ended without recording it (killed)
    return payload


def find_live_round(run: Path) -> tuple[Path, dict[str, object]]:
    """The live status of the running round, or of the newest round with snapshots."""
    rounds = []
    for path in (run / "AutoRefine").glob(f"round_*/{LIVE_DIRECTORY}/{STATUS_FILE}"):
        match = re.fullmatch(r"round_(\d+)", path.parent.parent.name)
        if match:
            rounds.append((int(match.group(1)), path))
    statuses = []
    for _, path in sorted(rounds, reverse=True):
        try:
            statuses.append((path, read_live_status(path)))
        except LiveViewError:
            continue
    for path, status in statuses:
        if status.get("state") == "running":
            return path, status
    for path, status in statuses:
        if status.get("snapshots"):
            return path, status
    if statuses:
        path, status = statuses[0]
        raise LiveViewError(
            f"{status.get('checkpoint')} has no live snapshot (state: {status.get('state')})"
        )
    raise LiveViewError(
        "No AutoRefine round in this run has live snapshots. Rounds started with this "
        "NASolve version save one when phenix.refine starts refining and after every "
        "macro-cycle."
    )


def latest_snapshot(status_path: Path, status: Mapping[str, object]) -> dict[str, object]:
    snapshots = status.get("snapshots")
    if not isinstance(snapshots, list) or not snapshots:
        raise LiveViewError(
            f"{status.get('checkpoint')} is {status.get('state')} but has no snapshot yet; "
            "phenix.refine saves the first one when it starts refining. Try again in a few seconds."
        )
    latest = snapshots[-1]
    if not isinstance(latest, dict):
        raise LiveViewError(f"{status_path} has a malformed snapshot")
    directory = status_path.parent
    model = directory / _plain_name(latest.get("model"), "snapshot model")
    maps = directory / _plain_name(latest.get("maps"), "snapshot maps")
    for path in (model, maps):
        if not path.is_file():
            raise LiveViewError(f"Snapshot file is missing: {path}")
    return {**latest, "model_path": model, "maps_path": maps}


def resolve_listed_path(value: object, live_directory: Path) -> Path | None:
    """A path written by ``LiveSnapshots`` (relative to the live directory, or absolute)."""
    if not isinstance(value, str) or not value:
        return None
    path = Path(value)
    return (path if path.is_absolute() else live_directory / path).resolve()


__all__ = [
    "GUI_ENVIRONMENT_VARIABLE",
    "LIVE_DIRECTORY",
    "LiveSnapshots",
    "LiveViewError",
    "MAP_LABELS",
    "PHENIX_STATUS_FILE",
    "STATUS_FILE",
    "find_live_round",
    "format_snapshot",
    "latest_snapshot",
    "read_live_status",
    "resolve_listed_path",
    "run_with_snapshots",
    "snapshot_label",
]
