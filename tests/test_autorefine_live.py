"""Watching a running AutoRefine round in Coot (user's request, 2026-11-26).

phenix.refine (Phenix 2.2.1) saves ``.status.pkl`` (model plus 2mFo-DFc and
mFo-DFc map coefficients) at the start, after every macro-cycle and at the
end when ``PHENIX_GUI_ENVIRONMENT`` is set: the hook its GUI uses. AutoRefine
sets it, converts each snapshot with phenix.python into
``round_NNN/live/snapshot_NNN.pdb`` and ``_maps.mtz``, lists them in
``live/status.json`` and prints R-work/R-free per snapshot; ``nasolve show RUN
--live`` opens the newest one in Coot with a button that loads the next.

phenix.refine, phenix.python and Coot are simulated: the stand-in refine writes
JSON "snapshots" and the stand-in phenix.python converts them.
"""

import json
import os
import stat
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import nasolve.live_refinement as live_refinement
from nasolve.autorefine import execute_autorefine
from nasolve.cli import build_parser
from nasolve.coot_view import CootViewError, _live_script, launch_coot_view

from .test_autorefine import PHENIX_21, make_mtz_dump, make_refine, make_refine_run


def snapshot_refine(root: Path, cycles: int, pause: float = 0.3) -> Path:
    """The stand-in phenix.refine, saving snapshots as phenix.refine does when
    PHENIX_GUI_ENVIRONMENT is set (start, each macro-cycle, end)."""
    inner = make_refine(root)
    wrapper = root / "phenix.refine.snapshots"
    wrapper.write_text(
        f"#!{sys.executable}\n"
        "import json, os, subprocess, sys, time\n"
        "from pathlib import Path\n"
        "gui = os.environ.get('PHENIX_GUI_ENVIRONMENT')\n"
        "Path('gui_env_dry_run.txt' if '--dry_run' in sys.argv else 'gui_env.txt')"
        ".write_text(repr(gui))\n"
        "def save(n, r_work):\n"
        "    if gui:\n"
        "        Path('.status.pkl').write_text(json.dumps({'n': n, 'r_work': r_work, "
        "'r_free': r_work + 0.03}))\n"
        f"    time.sleep({pause})\n"
        "if '--dry_run' not in sys.argv:\n"
        "    save(1, 0.40)\n"
        f"    for cycle in range(1, {cycles} + 1):\n"
        f"        print('REFINEMENT MACRO_CYCLE %d OF %d' % (cycle, {cycles}), flush=True)\n"
        "        save(cycle + 1, 0.40 - 0.05 * cycle)\n"
        f"    save({cycles} + 2, 0.40 - 0.05 * {cycles} - 0.001)\n"
        f"sys.exit(subprocess.call([{str(inner)!r}, *sys.argv[1:]]))\n"
    )
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)
    return wrapper


def fake_phenix_python(root: Path, *, fail: bool = False) -> Path:
    """Stand-in phenix.python: converts a JSON 'snapshot' when given the live
    converter; anything else goes to the stand-in make_refine wrote."""
    executable = root / "phenix.python"
    original = root / "phenix.python.orig"
    executable.rename(original)
    executable.write_text(
        f"#!{sys.executable}\n"
        "import json, subprocess, sys\n"
        "from pathlib import Path\n"
        "if len(sys.argv) < 3 or sys.argv[1] != '-c' or 'NASOLVE_LIVE_SNAPSHOT' not in sys.argv[2]:\n"
        f"    sys.exit(subprocess.call([{str(original)!r}, *sys.argv[1:]]))\n"
        + ("print('Traceback: cannot unpickle'); sys.exit(1)\n" if fail else "")
        + "source, model, maps = sys.argv[3:6]\n"
        "data = json.loads(Path(source).read_text())\n"
        "Path(model).write_text('REMARK snapshot %d\\nATOM      1  P    DC A   1       0.000   0.000   0.000  1.00 20.00           P\\n' % data['n'])\n"
        "Path(maps).write_bytes(b'maps %d' % data['n'])\n"
        "print('NASOLVE_LIVE_SNAPSHOT:' + json.dumps({'r_work': data['r_work'], 'r_free': data['r_free'], "
        "'bonds': 0.01, 'angles': 1.5}))\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IEXEC)
    return executable


class LiveRoundTests(unittest.TestCase):
    def setUp(self):
        patcher = mock.patch.object(live_refinement, "POLL_SECONDS", 0.05)
        patcher.start()
        self.addCleanup(patcher.stop)

    def refine(self, root: Path, run: Path, refine: Path, lines: list[str], cycles: int = 2):
        return execute_autorefine(
            run, refine, make_mtz_dump(root),
            phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
            macro_cycles=cycles, auto_select_success=False, live_report=lines.append,
        )

    def test_each_macro_cycle_is_saved_converted_and_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            refine = snapshot_refine(root, cycles=2)
            fake_phenix_python(root)
            lines: list[str] = []
            result = self.refine(root, run, refine, lines)
            self.assertEqual(result.status, "AUTOREFINE_READY", result.message)
            round_directory = result.report_path.parent
            # The hook is set for the refinement, not for the dry run.
            self.assertEqual((round_directory / "gui_env.txt").read_text(), "'1'")
            self.assertEqual((round_directory / "gui_env_dry_run.txt").read_text(), "None")
            status = json.loads((round_directory / "live" / "status.json").read_text())
            self.assertEqual(status["state"], "finished")
            self.assertEqual(status["checkpoint"], result.checkpoint_id)
            self.assertEqual(status["parent_checkpoint"], "postmr")
            self.assertEqual(status["macro_cycles"], 2)
            self.assertEqual(
                [(s["number"], s["label"], s["macro_cycle"]) for s in status["snapshots"]],
                [(1, "start", 0), (2, "macro-cycle 1 of 2", 1), (3, "macro-cycle 2 of 2", 2),
                 (4, "final", 2)])
            self.assertEqual([round(s["r_work"], 4) for s in status["snapshots"]],
                             [0.40, 0.35, 0.30, 0.299])
            for snapshot in status["snapshots"]:
                self.assertTrue((round_directory / "live" / snapshot["model"]).is_file())
                self.assertTrue((round_directory / "live" / snapshot["maps"]).is_file())
            self.assertEqual(sorted(p.name for p in (round_directory / "live").iterdir()), [
                "snapshot_001.pdb", "snapshot_001_maps.mtz", "snapshot_002.pdb",
                "snapshot_002_maps.mtz", "snapshot_003.pdb", "snapshot_003_maps.mtz",
                "snapshot_004.pdb", "snapshot_004_maps.mtz", "status.json"])
            # The start model and ligand dictionaries Coot needs, relative to live/.
            start = (round_directory / "live" / status["start_model"]).resolve()
            self.assertEqual(start, Path(json.loads(result.report_path.read_text())["inputs"]["model"]))
            self.assertTrue(all((round_directory / "live" / p).resolve().suffix == ".cif"
                                for p in status["dictionaries"]))
            self.assertEqual(lines, [
                f"  {result.checkpoint_id} start: R-work 0.4000, R-free 0.4300",
                f"  {result.checkpoint_id} macro-cycle 1 of 2: R-work 0.3500, R-free 0.3800",
                f"  {result.checkpoint_id} macro-cycle 2 of 2: R-work 0.3000, R-free 0.3300",
                f"  {result.checkpoint_id} final: R-work 0.2990, R-free 0.3290",
            ])
            record = json.loads(result.report_path.read_text())["live"]
            self.assertEqual(record["enabled"], True)
            self.assertEqual((record["state"], record["snapshots"], record["warnings"]),
                             ("finished", 4, []))
            self.assertEqual(record["status_file"], "live/status.json")

    def test_without_phenix_python_the_round_runs_as_before(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            refine = snapshot_refine(root, cycles=1, pause=0.0)
            (root / "phenix.python").unlink()
            lines: list[str] = []
            result = self.refine(root, run, refine, lines)
            self.assertEqual(result.status, "AUTOREFINE_READY", result.message)
            round_directory = result.report_path.parent
            self.assertEqual((round_directory / "gui_env.txt").read_text(), "None")
            self.assertFalse((round_directory / "live").exists())
            record = json.loads(result.report_path.read_text())["live"]
            self.assertEqual((record["enabled"], record["state"]), (False, "unavailable"))
            self.assertIn("phenix.python not found", record["reason"])
            self.assertEqual(lines, [])

    def test_failed_conversions_stop_the_live_view_but_not_the_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            refine = snapshot_refine(root, cycles=3)
            fake_phenix_python(root, fail=True)
            lines: list[str] = []
            result = self.refine(root, run, refine, lines)
            self.assertEqual(result.status, "AUTOREFINE_READY", result.message)
            record = json.loads(result.report_path.read_text())["live"]
            self.assertEqual((record["state"], record["snapshots"]), ("finished", 0))
            self.assertEqual(len(record["warnings"]), 4)
            self.assertIn("snapshot 1 (start) not converted: phenix.python exited 1: "
                          "Traceback: cannot unpickle", record["warnings"][0])
            self.assertIn("live view stopped after 3 failed conversions", record["warnings"][3])
            self.assertTrue(lines[0].startswith(f"  {result.checkpoint_id} live view: snapshot 1"))

    def test_ctrl_c_stops_phenix_refine_and_marks_the_live_view(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            refine = snapshot_refine(root, cycles=5)
            fake_phenix_python(root)
            lines: list[str] = []
            started: list[subprocess.Popen] = []
            real_popen = subprocess.Popen

            def popen(*args, **kwargs):
                process = real_popen(*args, **kwargs)
                started.append(process)
                return process

            calls = {"n": 0}
            real_poll = live_refinement.LiveSnapshots.poll

            def poll(self):
                calls["n"] += 1
                if calls["n"] == 30:
                    raise KeyboardInterrupt
                return real_poll(self)

            with mock.patch.object(live_refinement.subprocess, "Popen", popen), \
                    mock.patch.object(live_refinement.LiveSnapshots, "poll", poll):
                with self.assertRaises(KeyboardInterrupt):
                    self.refine(root, run, refine, lines)
            self.assertIsNotNone(started[-1].returncode)   # phenix.refine is not left running
            status_files = list((run / "AutoRefine").glob("round_*/live/status.json"))
            self.assertEqual(len(status_files), 1)
            self.assertEqual(json.loads(status_files[0].read_text())["state"], "interrupted")


class LiveViewTests(unittest.TestCase):
    def setUp(self):
        patcher = mock.patch.object(live_refinement, "POLL_SECONDS", 0.05)
        patcher.start()
        self.addCleanup(patcher.stop)

    def finished_round(self, root: Path) -> tuple[Path, object]:
        run = make_refine_run(root)
        refine = snapshot_refine(root, cycles=1, pause=0.2)
        fake_phenix_python(root)
        result = execute_autorefine(
            run, refine, make_mtz_dump(root),
            phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
            macro_cycles=1, auto_select_success=False,
        )
        self.assertEqual(result.status, "AUTOREFINE_READY", result.message)
        return run, result

    def launch(self, run: Path):
        captured = {}

        def launcher(command, **kwargs):
            captured["command"] = command
            captured["cwd"] = kwargs["cwd"]
            return SimpleNamespace(pid=4321)

        result = launch_coot_view(run, Path("/usr/bin/true"), live=True, launcher=launcher)
        return result, captured

    def test_live_view_opens_the_newest_snapshot_with_the_start_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, refined = self.finished_round(root)
            live = refined.report_path.parent / "live"
            result, captured = self.launch(run)
            command = captured["command"]
            self.assertEqual(command[1:3], ["--pdb", str((live / "snapshot_003.pdb").resolve())])
            start = Path(json.loads(refined.report_path.read_text())["inputs"]["model"])
            self.assertEqual(command[3:5], ["--pdb", str(start)])
            self.assertIn("--dictionary", command)
            self.assertEqual(command[command.index("--auto") + 1],
                             str((live / "snapshot_003_maps.mtz").resolve()))
            script = Path(command[command.index("--script") + 1]).read_text()
            self.assertIn("def nasolve_live_load_latest", script)
            self.assertNotIn("load_anomalous_map", script)
            self.assertEqual(result.stage, "autorefine-live")
            self.assertEqual(result.checkpoint_id, refined.checkpoint_id)
            self.assertIn("snapshot 3 of", result.source)
            self.assertIn("final; R-work 0.3490, R-free 0.3790; refinement finished", result.source)
            self.assertEqual(captured["cwd"], run / "CootGUI" / "autorefine-live" / refined.checkpoint_id)
            record = json.loads((captured["cwd"] / "launch.json").read_text())
            self.assertIsNone(record["anomalous_map"])
            self.assertEqual(record["live"]["number"], 3)

    def test_the_running_round_is_preferred_and_a_dead_one_is_not_running(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, refined = self.finished_round(root)
            # A later round still running (this process) with one snapshot.
            later = run / "AutoRefine" / "round_009" / "live"
            later.mkdir(parents=True)
            (later / "snapshot_001.pdb").write_text("REMARK running\n")
            (later / "snapshot_001_maps.mtz").write_bytes(b"maps")
            status = json.loads((refined.report_path.parent / "live" / "status.json").read_text())
            status.update(checkpoint="refine-009", state="running", pid=os.getpid(),
                          snapshots=[{"number": 1, "label": "start", "macro_cycle": 0,
                                      "model": "snapshot_001.pdb", "maps": "snapshot_001_maps.mtz",
                                      "r_work": 0.41, "r_free": 0.44}],
                          start_model="../../missing_start.pdb", dictionaries=[])
            (later / "status.json").write_text(json.dumps(status))
            result, captured = self.launch(run)
            self.assertEqual(result.checkpoint_id, "refine-009")
            self.assertIn("refinement running", result.source)
            self.assertEqual(result.extra_model_paths, ())   # start model missing: not opened
            # NASolve died without recording it: "stopped", still the newest with snapshots.
            dead = subprocess.Popen([sys.executable, "-c", "pass"])
            dead.wait()
            status["pid"] = dead.pid
            (later / "status.json").write_text(json.dumps(status))
            result, _ = self.launch(run)
            self.assertIn("refinement stopped", result.source)

    def test_nothing_to_show(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_refine_run(root)
            with self.assertRaisesRegex(CootViewError, "No AutoRefine round in this run has live"):
                self.launch(run)
            live = run / "AutoRefine" / "round_001" / "live"
            live.mkdir(parents=True)
            (live / "status.json").write_text(json.dumps({
                "schema_version": 1, "kind": "nasolve-autorefine-live", "checkpoint": "refine-001",
                "state": "running", "pid": os.getpid(), "snapshots": []}))
            with self.assertRaisesRegex(CootViewError, "refine-001 is running but has no "
                                                       "snapshot yet"):
                self.launch(run)
            with self.assertRaisesRegex(CootViewError, "cannot be combined"):
                launch_coot_view(run, Path("/usr/bin/true"), live=True, stage="autorefine",
                                 launcher=lambda *a, **k: SimpleNamespace(pid=1))

    def test_show_has_a_live_option(self):
        args = build_parser().parse_args(["show", "AutoMR/run_001", "--live"])
        self.assertTrue(args.live)
        self.assertFalse(build_parser().parse_args(["show"]).live)


class FakeCoot:
    """Enough of Coot's Python API for the live loader."""

    def __init__(self):
        self.molecules = {}   # imol -> dict(kind, path, diff, level, name)
        self.refinement_map = -1
        self.scroll_map = -1
        self.messages = []
        self.recentre = 1

    def add(self, kind, path, diff=False, level=None):
        imol = len(self.molecules)
        self.molecules[imol] = {"kind": kind, "path": path, "diff": diff, "level": level,
                                "name": path}
        return imol

    # API
    def graphics_n_molecules(self):
        return len(self.molecules)

    def is_valid_model_molecule(self, imol):
        return self.molecules.get(imol, {}).get("kind") == "model"

    def is_valid_map_molecule(self, imol):
        return self.molecules.get(imol, {}).get("kind") == "map"

    def molecule_name(self, imol):
        return self.molecules[imol]["path"]

    def mtz_hklin_for_map(self, imol):
        return self.molecules[imol]["path"]

    def set_molecule_name(self, imol, name):
        self.molecules[imol]["name"] = name

    def map_is_difference_map(self, imol):
        return 1 if self.molecules[imol]["diff"] else 0

    def get_contour_level_in_sigma(self, imol):
        return self.molecules[imol]["level"]

    def set_contour_level_in_sigma(self, imol, level):
        self.molecules[imol]["level"] = level

    def imol_refinement_map(self):
        return self.refinement_map

    def set_imol_refinement_map(self, imol):
        self.refinement_map = imol

    def scroll_wheel_map(self):
        return self.scroll_map

    def set_scroll_wheel_map(self, imol):
        self.scroll_map = imol

    def set_recentre_on_read_pdb(self, value):
        self.recentre = value

    def read_pdb(self, path):
        return self.add("model", path)

    def make_and_draw_map(self, path, amplitude, phase, weight, use_weights, is_diff):
        return self.add("map", path, diff=bool(is_diff), level=None)

    def close_molecule(self, imol):
        self.molecules[imol]["kind"] = "closed"

    def add_status_bar_text(self, text):
        self.messages.append(text)

    def info_dialog(self, text):
        self.messages.append("DIALOG " + text)


class LiveScriptTests(unittest.TestCase):
    def test_load_latest_replaces_the_snapshot_and_keeps_contours(self):
        with tempfile.TemporaryDirectory() as directory:
            live = Path(directory)
            status = {"state": "running", "snapshots": [
                {"number": 1, "label": "start", "model": "snapshot_001.pdb",
                 "maps": "snapshot_001_maps.mtz", "r_work": 0.41, "r_free": 0.44}]}
            (live / "status.json").write_text(json.dumps(status))
            coot = FakeCoot()
            first_model = coot.add("model", str(live / "snapshot_001.pdb"))
            start = coot.add("model", str(live / "start.pdb"))
            two = coot.add("map", str(live / "snapshot_001_maps.mtz"), level=1.4)
            diff = coot.add("map", str(live / "snapshot_001_maps.mtz"), diff=True, level=2.5)
            coot.refinement_map = two
            script = _live_script({
                "status": str(live / "status.json"), "checkpoint": "refine-009", "number": 1,
                "label": "start", "summary": "snapshot 1", "model": str(live / "snapshot_001.pdb"),
                "maps": str(live / "snapshot_001_maps.mtz"), "start_model": str(live / "start.pdb")})
            module = types.ModuleType("coot")
            for name in dir(coot):
                if not name.startswith("_"):
                    setattr(module, name, getattr(coot, name))
            namespace = {}
            with mock.patch.dict(sys.modules, {"coot": module, "gi": None, "coot_gui": None}):
                exec(script, namespace)
                self.assertEqual(coot.molecules[first_model]["name"], "Live snapshot 1 (start)")
                self.assertEqual(coot.molecules[start]["name"], "Start model of refine-009")
                self.assertEqual(coot.molecules[diff]["name"], "Live mFo-DFc (snapshot 1)")
                # Nothing newer yet.
                self.assertIsNone(namespace["nasolve_live_load_latest"]())
                self.assertIn("already showing the newest snapshot, 1 (start)", coot.messages[-1])
                status["snapshots"].append({"number": 2, "label": "macro-cycle 1 of 5",
                                            "model": "snapshot_002.pdb", "maps": "snapshot_002_maps.mtz",
                                            "r_work": 0.35, "r_free": 0.39})
                (live / "status.json").write_text(json.dumps(status))
                self.assertEqual(namespace["nasolve_live_load_latest"](), 2)
            new = {imol: m for imol, m in coot.molecules.items() if m["kind"] != "closed"}
            self.assertEqual(sorted(m["path"].rsplit("/", 1)[-1] for m in new.values()),
                             ["snapshot_002.pdb", "snapshot_002_maps.mtz", "snapshot_002_maps.mtz",
                              "start.pdb"])
            levels = {m["diff"]: m["level"] for m in new.values() if m["kind"] == "map"}
            self.assertEqual(levels, {False: 1.4, True: 2.5})   # the user's contour levels kept
            self.assertEqual(coot.molecules[coot.refinement_map]["path"],
                             str(live / "snapshot_002_maps.mtz"))
            self.assertEqual(coot.recentre, 0)   # the view stays where it is
            self.assertIn("snapshot 2 (macro-cycle 1 of 5): R-work 0.3500, R-free 0.3900; "
                          "refinement running", coot.messages[-1])


if __name__ == "__main__":
    unittest.main()
