"""AutoFDP step 2: one f''-only round from the step-1 model.

phenix.refine, phenix.mtz.dump, phenix.python and Coot are simulated. The
stand-in phenix.refine reads the anomalous groups NASolve wrote into its
parameters and echoes them, with a chosen refined f'', in phenix.refine's log
format, so the round, the f'' comparison and the record are exercised end to
end.
"""

import contextlib
import io
import json
import stat
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.autofdp import (
    STATUS_FAILED,
    STATUS_ROUND_DONE,
    AutoFDPError,
    fdp_round,
    held_parameter_changes,
    parse_theory_overrides,
    prepare_start_model,
    resolve_fdp_start,
)
from nasolve.checkpoints import list_checkpoints, read_registry
from nasolve.cli import main
from nasolve.config import AppConfig
from nasolve.coot_runtime import CootInstallation
from nasolve.phenix_runtime import PhenixInstallation

from .test_autofdp_bstart import make_refined_run
from .test_autorefine import PHENIX_21, make_mtz_dump

ENV = {"PATH": "/usr/bin:/bin"}
SILVER = "chain B and resid 101 and name AG"


def make_metal_run(root: Path) -> Path:
    """refine-001 holds the bridged Ag model (Ag B:101, occupancy 0.90,
    B 379.99, between A:12 DC with B 60 and B:4 DU with B 40); PostMR
    recorded the Ag as its one anomalous scatterer."""
    run = make_refined_run(root)
    report_path = run / "report.json"
    report = json.loads(report_path.read_text())
    report["postmr"]["mutation_actions"] = [
        {"site": "A:12", "before": "DT", "after": "DC"},
        {"site": "B:4", "before": "DT", "after": "DU"},
    ]
    report["postmr"]["anomalous"] = {
        "candidates": [],
        "metal_sites": [{
            "site": "B:101", "residue": "AG", "atom_name": "AG", "alternate": None,
            "element": "AG", "element_source": "element column",
            "coordinates": [0.0, 0.0, 0.0], "kind": "metal",
        }],
    }
    report_path.write_text(json.dumps(report))
    return run


def make_fdp_refine(root: Path, fdp: dict[str, float] | None = None, *,
                    change_b: bool = False) -> Path:
    """Stand-in phenix.refine that echoes each anomalous group it was given,
    with f'' from ``fdp`` (groups not in ``fdp`` are not reported)."""
    fdp = {SILVER: 2.6} if fdp is None else fdp
    executable = root / "phenix.refine"
    edit = (
        "lines = Path('refined_001.pdb').read_text().splitlines(keepends=True)\n"
        "lines[0] = lines[0][:60] + '%6.2f' % (float(lines[0][60:66]) + 5) + lines[0][66:]\n"
        "Path('refined_001.pdb').write_text(''.join(lines))\n"
        if change_b else ""
    )
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, re, shutil, sys\n"
        "from pathlib import Path\n"
        "if '--dry_run' in sys.argv:\n"
        "    print('Dry-run parameters accepted')\n"
        "    raise SystemExit(0)\n"
        "Path('received_args.json').write_text(json.dumps(sys.argv[1:]))\n"
        "params = next(Path(a) for a in sys.argv[1:] if a.endswith('.params')).read_text()\n"
        "model = next(Path(a) for a in sys.argv[1:] if a.endswith('.pdb'))\n"
        "shutil.copyfile(model, 'refined_001.pdb')\n"
        + edit +
        "Path('refined_001.cif').write_text('data_model\\n')\n"
        "Path('refined_001.reflections.cif').write_text('data_reflections\\n')\n"
        "Path('refined_001_map_coeffs.mtz').write_bytes(b'maps')\n"
        "Path('refined_001_data.mtz').write_bytes(b'refine evidence')\n"
        "print('Start: r_work = 0.2165 r_free = 0.2124 bonds = 0.012 angles = 1.82')\n"
        "print('Final: r_work = 0.2161 r_free = 0.2121 bonds = 0.012 angles = 1.82')\n"
        "print('clashscore = 4.8')\n"
        f"fdp = {fdp!r}\n"
        "for block in re.findall(r'group \\{(.*?)\\}', params, re.S):\n"
        "    selection = re.search(r'selection = \"(.*?)\"', block).group(1)\n"
        "    if selection not in fdp:\n"
        "        continue\n"
        "    print('Anomalous scatterer group:')\n"
        "    print(f'  Selection: \"{selection}\"')\n"
        "    print('  Number of selected scatterers: 1')\n"
        "    print('  f_prime: ' + re.search(r'f_prime = (\\S+)', block).group(1))\n"
        "    print(f'  f_double_prime: {fdp[selection]}')\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    python = root / "phenix.python"
    python.write_text(
        "#!/bin/sh\n"
        "case \"$*\" in\n"
        "  *NASOLVE_SCATTERING_JSON*) printf '%s\\n' 'NASOLVE_SCATTERING_JSON:{\"AG\": [-0.12, 2.03]}';;\n"
        "  *) printf '%s\\n' '2.03';;\n"
        "esac\n"
    )
    python.chmod(python.stat().st_mode | stat.S_IXUSR)
    return executable


def _round(run, root, *, start="autofdp-001", **keywords):
    return fdp_round(
        run, keywords.pop("refine", None) or make_fdp_refine(root),
        make_mtz_dump(root, wavelength=1.37762), phenix_version=PHENIX_21,
        start_checkpoint=start, environment=ENV, processor_count=1, **keywords,
    )


class FdpRoundTests(unittest.TestCase):
    def test_round_refines_only_the_anomalous_scatterers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            prepare_start_model(run, from_checkpoint="refine-001")
            result = _round(run, root)
            self.assertEqual(result.status, STATUS_ROUND_DONE, result.message)
            self.assertEqual(result.checkpoint_id, "refine-002")
            self.assertEqual(result.start_checkpoint, "autofdp-001")

            # Every strategy off except group_anomalous; f' fixed, f'' from 0.
            params = (result.round_directory / "autorefine.params").read_text()
            self.assertIn("strategy = *group_anomalous\n", params)
            self.assertIn(f'selection = "{SILVER}"', params)
            self.assertIn("f_prime = -0.12", params)
            self.assertIn("f_double_prime = 0.0", params)
            self.assertIn("refine = f_prime *f_double_prime", params)
            self.assertIn("optimize_xyz_weight = False", params)
            self.assertIn("optimize_adp_weight = False", params)
            self.assertIn("number_of_macro_cycles = 3", params)

            site = result.sites[0]
            self.assertTrue(site.autofdp_site)
            self.assertEqual((site.site, site.atom_name, site.element), ("B:101", "AG", "AG"))
            self.assertEqual(site.occupancy, 0.9)
            self.assertEqual(site.b_factor, 50.0)  # set by step 1
            self.assertEqual(site.f_prime, -0.12)
            self.assertEqual(site.f_double_prime_refined, 2.6)
            self.assertEqual(site.f_double_prime_theory, 2.03)
            self.assertEqual(site.theory_source, "cctbx-henke")
            self.assertAlmostEqual(site.deviation, 0.57)
            self.assertTrue(site.within_tolerance)
            self.assertEqual(result.warnings, ())

            # A normal AutoRefine checkpoint; the current one is unchanged.
            run_dir, registry = read_registry(run)
            self.assertEqual(registry["current"], "postmr")
            item = next(c for c in registry["checkpoints"] if c["id"] == "refine-002")
            self.assertEqual((item["parent"], item["recipe"]), ("autofdp-001", "AutoFDP/fdp"))

            record = json.loads(result.report_path.read_text())
            self.assertEqual(record["status"], STATUS_ROUND_DONE)
            entry = record["rounds"][0]
            self.assertEqual((entry["round"], entry["checkpoint"]), (1, "refine-002"))
            self.assertEqual(entry["strategy"], ["group_anomalous"])
            self.assertEqual(entry["starting_groups"][SILVER], {"f_prime": -0.12, "f_double_prime": 0.0})
            self.assertEqual(entry["sites"][0]["deviation"], 0.57)
            self.assertEqual(entry["held_parameter_changes"]["max_b_change"], 0.0)
            self.assertEqual(entry["wavelength"], 1.37762)
            log = (result.record_directory / "autofdp.log").read_text()
            self.assertIn("AutoFDP step 1", log)
            self.assertIn("AutoFDP step 2", log)
            self.assertIn("deviation +0.57", log)

            # A second round from the same B start adds to the same record.
            again = _round(run, root, refine=make_fdp_refine(root, {SILVER: 5.0}))
            self.assertEqual(again.checkpoint_id, "refine-003")
            self.assertEqual(again.record_directory, result.record_directory)
            self.assertAlmostEqual(again.sites[0].deviation, 2.97)
            self.assertFalse(again.sites[0].within_tolerance)
            self.assertIn("1 of 1 site(s) more than 1 e- from theory", again.message)
            self.assertEqual(len(json.loads(again.report_path.read_text())["rounds"]), 2)

    def test_theory_override_and_tolerance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            prepare_start_model(run, from_checkpoint="refine-001")
            result = _round(run, root, refine=make_fdp_refine(root, {SILVER: 5.0}),
                            theory_overrides={"AG": 4.5}, tolerance=0.4, cycles=2)
            site = result.sites[0]
            self.assertEqual(site.f_double_prime_theory, 4.5)
            self.assertEqual(site.theory_source, "user (--theory)")
            self.assertAlmostEqual(site.deviation, 0.5)
            self.assertFalse(site.within_tolerance)
            self.assertIn("number_of_macro_cycles = 2",
                          (result.round_directory / "autorefine.params").read_text())
            record = json.loads(result.report_path.read_text())
            self.assertEqual(record["rounds"][0]["theory"][SILVER]["f_double_prime_henke"], 2.03)
            # f' stays the calculated value even when f'' theory is overridden.
            self.assertEqual(record["rounds"][0]["starting_groups"][SILVER]["f_prime"], -0.12)
            with self.assertRaisesRegex(AutoFDPError, "HG, which is not an anomalous scatterer"):
                _round(run, root, theory_overrides={"HG": 10.17})

    def test_theory_override_parsing(self):
        self.assertEqual(parse_theory_overrides(["hg=10.17", "AG=2.06"]), {"HG": 10.17, "AG": 2.06})
        for bad in (["HG"], ["HG=x"], ["HG=-1"], ["H1=2"], ["HG=1", "HG=2"]):
            with self.subTest(bad=bad), self.assertRaises(AutoFDPError):
                parse_theory_overrides(bad)

    def test_starts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            self.assertTrue(resolve_fdp_start(run, "refine-001").needs_b_start)
            with self.assertRaisesRegex(AutoFDPError, "not an AutoRefine result or an AutoFDP B start"):
                resolve_fdp_start(run, "postmr")
            with self.assertRaisesRegex(AutoFDPError, "starts from the AutoFDP B start"):
                _round(run, root, start="refine-001")
            prepare_start_model(run, from_checkpoint="refine-001")
            start = resolve_fdp_start(run, "autofdp-001")
            self.assertFalse(start.needs_b_start)
            self.assertEqual(start.record_directory, run / "AutoFDP" / "fdp_001")

    def test_site_must_be_a_postmr_anomalous_scatterer(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            prepare_start_model(run, from_checkpoint="refine-001")
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["postmr"]["anomalous"]["metal_sites"] = []
            report["postmr"]["anomalous"]["candidates"] = [{
                "site": "B:4", "residue": "DU", "atom_name": "O4", "alternate": None,
                "element": "O", "coordinates": [-2.2, 0.0, 0.0],
            }]
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(AutoFDPError, "Not among the run's anomalous scatterers"):
                _round(run, root)
            record = json.loads((run / "AutoFDP" / "fdp_001" / "report.json").read_text())
            self.assertNotIn("rounds", record)
            self.assertEqual(len(list_checkpoints(run)), 3)  # no round was started

    def test_missing_refined_fdp_or_changed_held_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            prepare_start_model(run, from_checkpoint="refine-001")
            result = _round(run, root, refine=make_fdp_refine(root, {}))
            self.assertEqual(result.status, STATUS_FAILED)
            self.assertIn("no refined f''", result.message)
            self.assertIsNone(result.sites[0].deviation)
            record = json.loads(result.report_path.read_text())
            self.assertEqual(record["rounds"][0]["status"], STATUS_FAILED)

            changed = _round(run, root, refine=make_fdp_refine(root, change_b=True))
            self.assertIn("Held parameters changed", changed.warnings[0])
            self.assertEqual(
                held_parameter_changes(
                    run / "AutoFDP" / "fdp_001" / "start_model.pdb", changed.model_path
                )["max_b_change"],
                5.0,
            )


class AutoFDPCommandTests(unittest.TestCase):
    def _main(self, root: Path, arguments: list[str], launch=None):
        refine = make_fdp_refine(root)
        installation = PhenixInstallation(
            root=None, setup_script=None, version=PHENIX_21,
            executables={"phenix.refine": refine,
                         "phenix.mtz.dump": make_mtz_dump(root, wavelength=1.37762)},
            environment=ENV, source="test",
        )
        coot = root / "coot"
        coot.write_text("#!/bin/sh\n")
        launch = launch or (lambda *args, **kwargs: SimpleNamespace(
            pid=4321, model_path=Path("refined_001.pdb"), map_path=Path("refined_001_map_coeffs.mtz"),
            log_path=Path("coot_gui.log"),
        ))
        out, err = io.StringIO(), io.StringIO()
        with (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli.discover_phenix", return_value=installation),
            patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
            patch("nasolve.cli.discover_coot", return_value=CootInstallation(coot, "test", "fixture")),
            patch("nasolve.cli.remember_coot"),
            patch("nasolve.cli.launch_coot_view", side_effect=launch) as launched,
            contextlib.redirect_stdout(out), contextlib.redirect_stderr(err),
        ):
            code = main(["autofdp", *arguments])
        return code, out.getvalue(), err.getvalue(), launched

    def test_command_runs_both_steps_and_opens_coot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            code, out, err, launched = self._main(root, [str(run), "--from", "refine-001"])
            self.assertEqual(code, 0, err)
            self.assertIn("B:101 AG (AG), occupancy 0.90: B 379.99 -> 50.00", out)
            self.assertIn("Checkpoint: autofdp-001 (current checkpoint unchanged)", out)
            self.assertIn("f''-only round refine-002 from autofdp-001 (3 cycles", out)
            row = next(line for line in out.splitlines() if line.startswith("B:101 AG"))
            for value in ("0.90", "50.00", "-0.12", "2.60", "2.03", "+0.57", "within +/-1"):
                self.assertIn(value, row)
            self.assertIn("Theory: AG cctbx-henke", out)
            self.assertIn("AUTOFDP_FDP_ROUND_DONE", out)
            self.assertIn("Opened Coot (PID 4321) with refine-002", out)
            self.assertEqual(launched.call_args.kwargs["checkpoint"], "refine-002")
            self.assertIn("Occupancy adjustment and the later AutoFDP steps are not implemented yet", out)

            # Again from the B start, without Coot: a new round, same record.
            code, out, err, launched = self._main(
                root, [str(run), "--from", "autofdp-001", "--no-coot", "--theory", "AG=2.06"]
            )
            self.assertEqual(code, 0, err)
            launched.assert_not_called()
            self.assertIn("--checkpoint refine-003", out)
            self.assertIn("Theory: AG user (--theory)", out)
            self.assertNotIn("AutoFDP step 1", out)
            self.assertEqual(sorted(p.name for p in (run / "AutoFDP").iterdir()), ["fdp_001"])

    def test_coot_problem_does_not_fail_the_run(self):
        from nasolve.coot_view import CootViewError

        def broken(*args, **kwargs):
            raise CootViewError("fixture: no display")

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            code, out, err, _ = self._main(root, [str(run), "--from", "refine-001"], launch=broken)
            self.assertEqual(code, 0)
            self.assertIn("Coot view error: fixture: no display", err)
            self.assertIn("nasolve show", out)

    def test_plan_runs_no_phenix_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            out, err = io.StringIO(), io.StringIO()
            with (
                patch("nasolve.cli.load_config", return_value=AppConfig()),
                patch("nasolve.cli.discover_phenix", side_effect=AssertionError("no Phenix in --plan")),
                contextlib.redirect_stdout(out), contextlib.redirect_stderr(err),
            ):
                code = main(["autofdp", str(run), "--from", "refine-001", "--plan"])
            self.assertEqual(code, 0, err.getvalue())
            self.assertIn("B 379.99 -> 50.00", out.getvalue())
            self.assertIn("Then: one f''-only phenix.refine round of 3 cycles", out.getvalue())
            self.assertFalse((run / "AutoFDP").exists())

    def test_command_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_metal_run(root)
            code, _, err, _ = self._main(root, [str(run), "--from", "refine-001", "--theory", "AG"])
            self.assertEqual(code, 1)
            self.assertIn("--theory must look like HG=10.17", err)
            code, _, err, _ = self._main(root, [str(run), "--from", "refine-001", "--cycles", "0"])
            self.assertEqual(code, 1)
            self.assertFalse((run / "AutoFDP").exists())
            self._main(root, [str(run), "--from", "refine-001", "--no-coot"])
            code, _, err, _ = self._main(root, [str(run), "--from", "autofdp-001", "--site", "B:101"])
            self.assertEqual(code, 1)
            self.assertIn("already has its metal B set", err)

            # The round fails after step 1: the B start is kept and named.
            code, _, err, _ = self._main(
                root, [str(run), "--from", "refine-001", "--no-coot", "--theory", "HG=10.17"]
            )
            self.assertEqual(code, 1)
            self.assertIn("HG, which is not an anomalous scatterer", err)
            self.assertIn("The B start autofdp-002 (step 1) was kept", err)
            self.assertIn("--from autofdp-002", err)


if __name__ == "__main__":
    unittest.main()
