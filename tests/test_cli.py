import json
import io
import shlex
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.autorefine import AutoRefineError, AutoRefineResult
from nasolve.cli import _workspace_run, build_parser, main
from nasolve.config import AppConfig, ConfigError, WorkspaceSettings
from nasolve.coot_runtime import CootDiscoveryError, CootInstallation
from nasolve.coot_view import CootViewError, launch_coot_view
from nasolve.model_assessment import file_sha256
from nasolve.refine_doctor import RefineDoctorError

from .helpers import pdb_record, symlinked_temporary_directory
from .test_coot_view import add_manual_checkpoint, add_review_checkpoint, make_view_run
from .test_checkpoints import make_checkpoint_run


class CLITests(unittest.TestCase):
    def test_checkpoint_listing_cli_is_read_only_for_explicit_and_active_run(self):
        with tempfile.TemporaryDirectory() as directory:
            run = make_checkpoint_run(Path(directory))
            config = AppConfig(workspace=WorkspaceSettings(run=str(run)))
            with patch('nasolve.cli.load_config', return_value=config), patch('nasolve.cli.save_config') as save:
                for arguments in (['checkpoints', 'list', str(run)], ['checkpoints', 'list']):
                    with self.subTest(arguments=arguments), redirect_stdout(io.StringIO()) as output:
                        self.assertEqual(main(arguments), 0)
                        self.assertIn('postmr', output.getvalue())
                        self.assertFalse((run / 'AutoRefine').exists())
                save.assert_not_called()

    def test_autorefine_prints_exact_inspection_command_without_selecting(self):
        with tempfile.TemporaryDirectory(prefix='nasolve inspection ') as directory:
            run = make_checkpoint_run(Path(directory))
            report_path = run / 'refinement-report.json'
            report_path.write_text('{}')
            installation = SimpleNamespace(
                version='2.2.1', environment={},
                executables={'phenix.refine': Path('/phenix/refine'), 'phenix.mtz.dump': Path('/phenix/dump')},
            )
            for status, selected, products, viewable in (
                ('AUTOREFINE_READY', True, True, True),
                ('AUTOREFINE_REVIEW', False, True, True),
                ('AUTOREFINE_ANOMALOUS_FALLBACK', False, True, True),
                ('AUTOREFINE_FAILED', False, True, True),
                ('AUTOREFINE_FAILED', False, False, False),
                ('AUTOREFINE_REVIEW', False, False, False),
                # The resolver can find recorded maps absent from the result field.
                ('AUTOREFINE_REVIEW', False, False, True),
                # Non-null result paths do not establish artifact integrity.
                ('AUTOREFINE_REVIEW', False, True, False),
            ):
                result = AutoRefineResult(
                    status=status, message='fixture', exit_code=2 if status in {'AUTOREFINE_FAILED', 'AUTOREFINE_ANOMALOUS_FALLBACK'} else 0,
                    run_directory=run, round_directory=run / 'AutoRefine' / 'round_002',
                    checkpoint_id='refine-002', parent_checkpoint='postmr',
                    report_path=report_path, log_path=run / 'refine.log',
                    model_path=run / 'refined.pdb',
                    model_cif=None, reflection_cif=None,
                    map_coefficients=run / 'maps.mtz' if products else None,
                    statistics={}, selected_as_current=selected,
                )
                with (
                    self.subTest(status=status, products=products, viewable=viewable),
                    patch('nasolve.cli.load_config', return_value=AppConfig()),
                    patch('nasolve.cli.discover_phenix', return_value=installation),
                    patch('nasolve.cli.remember_phenix'), patch('nasolve.cli.save_config'),
                    patch('nasolve.cli.execute_autorefine', return_value=result),
                    patch(
                        'nasolve.cli.resolve_view_profile',
                        side_effect=None if viewable else CootViewError('fixture map unavailable'),
                    ) as resolve,
                    patch('nasolve.cli.select_checkpoint') as select,
                    redirect_stdout(io.StringIO()) as output,
                ):
                    self.assertEqual(main(['autorefine', str(run)]), result.exit_code)
                    text = output.getvalue()
                    commands = [line.removeprefix('Inspect: ') for line in text.splitlines() if line.startswith('Inspect: ')]
                    resolve.assert_called_once_with(run, checkpoint='refine-002')
                    if viewable:
                        self.assertEqual(len(commands), 1)
                        self.assertEqual(shlex.split(commands[0]), ['./nasolve', 'show', str(run), '--checkpoint', 'refine-002'])
                    else:
                        self.assertEqual(commands, [])
                        self.assertIn('Model/map inspection unavailable', text)
                        self.assertIn('fixture map unavailable', text)
                    self.assertIn('Current checkpoint updated.' if selected else 'Current checkpoint unchanged.', text)
                    self.assertNotIn('select this result', text)
                    select.assert_not_called()

    def test_show_active_run_and_explicit_overrides_report_exact_sources(self):
        with symlinked_temporary_directory() as directory:
            root = Path(directory)
            active = make_view_run(root, 1)
            newest = make_view_run(root, 2)
            add_review_checkpoint(active)
            config = AppConfig(workspace=WorkspaceSettings(
                dataset=str(active.parent.parent), run=str(active),
            ))
            coot = root / "coot"
            coot.write_text("")
            launch = partial(
                launch_coot_view,
                launcher=lambda *args, **kwargs: SimpleNamespace(pid=1234),
            )
            cases = (
                (["show"], active, "refine-001", "current refinement checkpoint refine-001"),
                (["show", str(active), "--checkpoint", "refine-002"], active,
                 "refine-002", "requested refinement checkpoint refine-002 (REVIEW)"),
                (["show", str(newest), "--stage", "autosol"], newest,
                 None, "PostMR ReadySet prepared model with separate AutoSol HA sites"),
                (["show", "last"], newest, "refine-001", "current refinement checkpoint refine-001"),
            )
            with (
                patch("nasolve.cli.load_config", return_value=config),
                patch("nasolve.cli.save_config"),
                patch("nasolve.cli.discover_coot", return_value=CootInstallation(coot, "test", "fixture")),
                patch("nasolve.cli.launch_coot_view", side_effect=launch),
            ):
                for arguments, run, checkpoint, source in cases:
                    output = io.StringIO()
                    with self.subTest(arguments=arguments), redirect_stdout(output):
                        self.assertEqual(main(arguments), 0)
                    text = output.getvalue()
                    self.assertIn(f"in run: {run.resolve()}", text)
                    self.assertIn(source, text)
                    self.assertIn("Map source:", text)
                    pen = run / "CootGUI" / ("autorefine" if checkpoint else "autosol")
                    if checkpoint:
                        pen /= checkpoint
                    record = json.loads((pen / "launch.json").read_text())
                    self.assertEqual(record["run_directory"], str(run.resolve()))
                    self.assertEqual(record["checkpoint_id"], checkpoint)
                    self.assertIn(record["map_source"], text)
                    self.assertIn(record["model_path"], text)
                    self.assertIn(record["map_path"], text)
                    self.assertEqual(config.workspace.run, str(active))

    def test_bare_show_opens_current_imported_model_in_active_workspace(self):
        with symlinked_temporary_directory() as directory:
            root = Path(directory)
            run = make_view_run(root)
            model = add_manual_checkpoint(run)
            config = AppConfig(workspace=WorkspaceSettings(run=str(run)))
            coot = root / "coot"
            coot.write_text("")
            output = io.StringIO()
            with (
                patch("nasolve.cli.load_config", return_value=config),
                patch("nasolve.cli.save_config"),
                patch("nasolve.cli.discover_coot", return_value=CootInstallation(coot, "test", "fixture")),
                patch("nasolve.cli.launch_coot_view", side_effect=partial(
                    launch_coot_view, launcher=lambda *args, **kwargs: SimpleNamespace(pid=1234),
                )),
                redirect_stdout(output),
            ):
                self.assertEqual(main(["show"]), 0)
            self.assertIn(f"Model: {model.resolve()}", output.getvalue())
            self.assertIn("current manual checkpoint manual-001", output.getvalue())
            self.assertIn("inherited refinement checkpoint refine-001", output.getvalue())

    def test_autorefine_forwards_discovered_phenix_version(self):
        installation = SimpleNamespace(
            version="1.20.1-4487",
            executables={
                "phenix.refine": Path("/phenix/phenix.refine"),
                "phenix.mtz.dump": Path("/phenix/phenix.mtz.dump"),
            },
            environment={"PATH": "/phenix"},
        )
        with (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli._workspace_run", return_value=Path("/dataset/AutoMR/run_001")),
            patch("nasolve.cli.discover_phenix", return_value=installation),
            patch("nasolve.cli.remember_phenix"),
            patch("nasolve.cli.save_config"),
            patch(
                "nasolve.cli.execute_autorefine",
                side_effect=AutoRefineError("sentinel"),
            ) as execute,
        ):
            self.assertEqual(main(["autorefine", "run_001"]), 2)

        self.assertEqual(execute.call_args.kwargs["phenix_version"], installation.version)

    def test_refine_doctor_forwards_discovered_phenix_version(self):
        installation = SimpleNamespace(
            version="2.1-6048",
            executables={
                "phenix.refine": Path("/phenix/phenix.refine"),
                "phenix.mtz.dump": Path("/phenix/phenix.mtz.dump"),
            },
            environment={"PATH": "/phenix"},
        )
        with (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli._workspace_run", return_value=Path("/dataset/AutoMR/run_001")),
            patch("nasolve.cli.discover_phenix", return_value=installation),
            patch("nasolve.cli.remember_phenix"),
            patch("nasolve.cli.save_config"),
            patch(
                "nasolve.cli.execute_refine_doctor",
                side_effect=RefineDoctorError("sentinel"),
            ) as execute,
        ):
            self.assertEqual(main(["refine-doctor", "run_001"]), 2)

        self.assertEqual(execute.call_args.kwargs["phenix_version"], installation.version)

    def test_refine_doctor_unknown_version_is_a_guarded_cli_error(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "dataset" / "AutoMR" / "run_001"
            installation = SimpleNamespace(
                version="unknown",
                executables={
                    "phenix.refine": Path("/phenix/phenix.refine"),
                    "phenix.mtz.dump": Path("/phenix/phenix.mtz.dump"),
                },
                environment={"PATH": "/phenix"},
            )
            with (
                patch("nasolve.cli.load_config", return_value=AppConfig()),
                patch("nasolve.cli._workspace_run", return_value=run),
                patch("nasolve.cli.discover_phenix", return_value=installation),
                patch("nasolve.cli.remember_phenix"),
                patch("nasolve.cli.save_config"),
            ):
                self.assertEqual(main(["refine-doctor", "run_001"]), 2)

            self.assertFalse((run / "RefineDoctor").exists())

    def test_workspace_drives_checkpoint_commands_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = root / "dataset"
            run = dataset / "AutoMR" / "run_004"
            model = run / "PostMR" / "Model" / "readyset_model.pdb"
            model.parent.mkdir(parents=True)
            model.write_text(pdb_record("ATOM", 1, "P", "DA", "A", 1) + "END\n")
            observations = dataset / "staraniso.mtz"
            observations.write_bytes(b"observations")
            report = {
                "workflow": "automr",
                "inputs": {
                    "reflections": str(observations),
                    "reflections_sha256": file_sha256(observations),
                },
                "postmr": {
                    "status": "POSTMR_READY",
                    "prepared_model": str(model),
                    "prepared_sha256": file_sha256(model),
                    "restraints": [],
                },
            }
            (run / "report.json").write_text(json.dumps(report))
            config_path = root / "config.json"

            with patch.dict(
                "os.environ", {"NASOLVE_CONFIG_FILE": str(config_path)}, clear=False
            ):
                self.assertEqual(main(["workspace", "use", str(run)]), 0)
                self.assertEqual(main(["checkpoints", "list"]), 0)
                self.assertEqual(main(["checkpoints", "use", "postmr"]), 0)
                self.assertEqual(main(["workspace", "clear"]), 0)

    def test_workspace_run_rejects_non_object_report(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "dataset" / "AutoMR" / "run_001"
            run.mkdir(parents=True)
            (run / "report.json").write_text(json.dumps([]))
            config = AppConfig(workspace=WorkspaceSettings(run=str(run)))

            with self.assertRaisesRegex(ConfigError, "not a JSON object"):
                _workspace_run(None, config)

    def test_workspace_and_optional_active_targets(self):
        parser = build_parser()
        workspace = parser.parse_args(["workspace", "use", "dataset/AutoMR/run_004"])
        self.assertEqual(workspace.workspace_action, "use")
        self.assertEqual(str(workspace.target), "dataset/AutoMR/run_004")
        self.assertIsNone(parser.parse_args(["autorefine"]).run)
        self.assertIsNone(parser.parse_args(["checkpoints", "list"]).run)
        self.assertIsNone(parser.parse_args(["show"]).target)

    def test_standard_frame_shortcuts(self):
        parser = build_parser()
        w = parser.parse_args(["automr", "dataset", "-W", "--pair", "D:T"])
        self.assertEqual((w.command, w.frame, w.pair), ("automr", "W", "D:T"))
        three_gbi = parser.parse_args(["automr", "dataset", "-3GBI", "--pair", "D:T"])
        self.assertEqual(three_gbi.frame, "3GBI")
        shunted = parser.parse_args([
            "automr", "dataset", "-W", "--pair", "D:T", "--allow-p1-standard",
            "--execute",
        ])
        self.assertTrue(shunted.allow_p1_standard)
        self.assertTrue(shunted.execute)
        mirrored = parser.parse_args([
            "automr", "dataset", "-W", "--pair", "D:T", "--mirror",
        ])
        self.assertTrue(mirrored.mirror)
        forced = parser.parse_args([
            "automr", "dataset", "-W", "--pair", "D:T",
            "--model", "5W6W_noPO4.pdb",
        ])
        self.assertEqual(forced.model, "5W6W_noPO4.pdb")

    def test_postmr_command(self):
        parser = build_parser()
        args = parser.parse_args([
            "postmr", "run_004", "--allow-mr-review", "--modified-pairs-only",
        ])
        self.assertEqual(args.command, "postmr")
        self.assertTrue(args.allow_mr_review)
        self.assertTrue(args.modified_pairs_only)

    def test_tfz_override_is_forwarded_and_reported(self):
        installation = SimpleNamespace(
            version="2.2.1",
            executables={
                "phenix.ready_set": Path("/phenix/phenix.ready_set"),
                "phenix.autosol": Path("/phenix/phenix.autosol"),
                "phenix.mtz.dump": Path("/phenix/phenix.mtz.dump"),
            },
            environment={"PATH": "/phenix"},
        )
        record = {"accepted_by": "--tfz-override", "mr_status": "MR_FAILED", "best_tfz": 5.1}
        warning = "MR_FAILED (best TFZ 5.10, below the 7.0 floor) because of --tfz-override"
        run = Path("/dataset/AutoMR/run_001")
        common = (
            patch("nasolve.cli.load_config", return_value=AppConfig()),
            patch("nasolve.cli._workspace_run", return_value=run),
            patch("nasolve.cli.discover_phenix", return_value=installation),
            patch("nasolve.cli.remember_phenix"),
            patch("nasolve.cli.save_config"),
            patch("nasolve.cli.discover_coot", side_effect=CootDiscoveryError("no Coot")),
        )
        postmr_result = SimpleNamespace(
            status="POSTMR_READY", message="ready", run_directory=run, input_model=None,
            model_path=run / "PostMR" / "model.pdb", metal_sites=(), anomalous_scatterer_phil=None,
            readyset_log=run / "PostMR" / "readyset.log", restraint_paths=(),
            report_path=run / "PostMR" / "report.json", tfz_override_record=record,
            narestraints_summary="NARestraints: 19 base pairs (guessed from the model) and stacking",
            overwritten=(),
        )
        metal_result = SimpleNamespace(
            status="AUTOSOL_METAL_NOT_FOUND", message="no sites", metals=("Ag",),
            partial_model_rmsd=1.0, report_path=Path("/missing/report.json"), mutations=(),
            partial_model=None, phaser_model=run / "Phaser" / "mr_solution.pdb",
            run_directory=run, output_directory=Path("/dataset/AutoSolMetal/run_001"),
            log_path=Path("/dataset/AutoSolMetal/run_001/autosol.log"), sites=(),
            tfz_override_record=record,
        )
        for arguments, target, result in (
            (["postmr", "run_001", "--tfz-override"], "nasolve.cli.prepare_postmr", postmr_result),
            (["autosol", "run_001", "--metal", "AG", "--tfz-override"],
             "nasolve.cli.execute_metal_search", metal_result),
        ):
            with self.subTest(arguments=arguments):
                output = io.StringIO()
                with (common[0], common[1], common[2], common[3], common[4], common[5],
                      patch(target, return_value=result) as call, redirect_stdout(output)):
                    self.assertEqual(main(arguments), 0)
                self.assertIs(call.call_args.kwargs["tfz_override"], True)
                self.assertIn(warning, output.getvalue())
                # Without the flag nothing is overridden.
                with (common[0], common[1], common[2], common[3], common[4], common[5],
                      patch(target, return_value=result) as call, redirect_stdout(io.StringIO())):
                    main(arguments[:-1])
                self.assertIs(call.call_args.kwargs["tfz_override"], False)

        # Without --metal the flag is refused, like the other metal-search options.
        error = io.StringIO()
        with redirect_stderr(error):
            self.assertEqual(main(["autosol", "run_001", "--tfz-override"]), 2)
        self.assertIn("--tfz-override", error.getvalue())

    def test_postmr_merged_model_is_forwarded_and_reported(self):
        installation = SimpleNamespace(
            version="2.2.1", executables={"phenix.ready_set": Path("/phenix/phenix.ready_set")},
            environment={"PATH": "/phenix"},
        )
        run = Path("/dataset/AutoMR/run_004")
        result = SimpleNamespace(
            status="POSTMR_READY", message="ready", run_directory=run,
            input_model=run / "Phaser" / "mr_solution-merged.pdb", input_model_source="metal-merged",
            merge_selected_by="--merged-model", other_merges=("mr_solution-merged-2.pdb",),
            merge_warnings=("mr_solution-merged.pdb has changed since it was merged; using the file as it is now",),
            model_path=run / "PostMR" / "model.pdb", metal_sites=(), anomalous_scatterer_phil=None,
            readyset_log=run / "PostMR" / "readyset.log", restraint_paths=(),
            report_path=run / "PostMR" / "report.json", tfz_override_record=None,
            narestraints_summary="NARestraints: not run; no base-pair or stacking restraints",
            overwritten=(),
        )
        output = io.StringIO()
        with (patch("nasolve.cli.load_config", return_value=AppConfig()),
              patch("nasolve.cli._workspace_run", return_value=run),
              patch("nasolve.cli.discover_phenix", return_value=installation),
              patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
              patch("nasolve.cli.discover_coot", side_effect=CootDiscoveryError("no Coot")),
              patch("nasolve.cli.prepare_postmr", return_value=result) as call,
              redirect_stdout(output)):
            self.assertEqual(main(["postmr", str(run), "--merged-model", "mr_solution-merged.pdb"]), 0)
        self.assertEqual(call.call_args.kwargs["merged_model"], "mr_solution-merged.pdb")
        self.assertIs(call.call_args.kwargs["use_merged_model"], True)
        text = output.getvalue()
        self.assertIn(f"Starting model: {run / 'Phaser' / 'mr_solution-merged.pdb'} (Phaser model with "
                      "merged metal sites, chosen with --merged-model)", text)
        self.assertIn("Other merged models recorded (choose one with --merged-model NAME): "
                      "mr_solution-merged-2.pdb", text)
        self.assertIn("Warning: mr_solution-merged.pdb has changed since it was merged", text)

    def test_autorefine_forwards_and_prints_extra_restraints(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root / "run_001"
            round_directory = run / "AutoRefine" / "round_001"
            round_directory.mkdir(parents=True)
            copy = round_directory / "extra_restraints" / "5W6W_secondary_structure.eff"
            report = round_directory / "report.json"
            report.write_text(json.dumps({"inputs": {"observation_labels": ["F"], "user_restraints": [
                {"source": "/data/5W6W_secondary_structure.eff", "sha256": "0" * 64, "copy": str(copy)},
            ]}, "refinement": {}}))
            result = SimpleNamespace(
                status="AUTOREFINE_READY", message="done", checkpoint_id="refine-001",
                parent_checkpoint="postmr", report_path=report, statistics={}, model_path=None,
                model_cif=None, reflection_cif=None, map_coefficients=None,
                log_path=round_directory / "phenix.refine.log", selected_as_current=False,
                run_directory=run, exit_code=2,
            )
            phenix = SimpleNamespace(version="2.2.1", environment=None, executables={
                "phenix.refine": Path("/phenix/phenix.refine"),
                "phenix.mtz.dump": Path("/phenix/phenix.mtz.dump"),
            })
            output = io.StringIO()
            with (patch("nasolve.cli.load_config", return_value=AppConfig()),
                  patch("nasolve.cli._workspace_run", return_value=run),
                  patch("nasolve.cli.discover_phenix", return_value=phenix),
                  patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
                  patch("nasolve.cli.resolve_view_profile", side_effect=CootViewError("no view")),
                  patch("nasolve.cli.execute_autorefine", return_value=result) as call,
                  redirect_stdout(output)):
                main(["autorefine", "run_001", "--restraints", "ss.eff", "--restraints", "x.cif"])
            self.assertEqual(call.call_args.kwargs["user_restraints"], [Path("ss.eff"), Path("x.cif")])
            self.assertIn("Extra restraints: 5W6W_secondary_structure.eff "
                          "(copy of /data/5W6W_secondary_structure.eff)", output.getvalue())

    def test_autosol_command(self):
        parser = build_parser()
        args = parser.parse_args(["autosol", "run_001"])
        self.assertEqual(args.command, "autosol")
        self.assertEqual(str(args.run), "run_001")

    def test_autorefine_and_checkpoint_commands(self):
        parser = build_parser()
        refine = parser.parse_args([
            "autorefine", "run_001", "--from", "clean", "--cycles", "5",
        ])
        self.assertEqual(refine.command, "autorefine")
        self.assertEqual(refine.from_checkpoint, "clean")
        self.assertEqual(refine.cycles, 5)
        self.assertEqual(refine.restraints, [])
        extra = parser.parse_args([
            "autorefine", "run_001", "--restraints", "a.eff", "--restraints", "DYI.cif",
        ])
        self.assertEqual(extra.restraints, [Path("a.eff"), Path("DYI.cif")])
        doctor = parser.parse_args([
            "refine-doctor", "run_001", "--from", "refine-003", "--cycles", "3",
        ])
        self.assertEqual(doctor.command, "refine-doctor")
        self.assertEqual(doctor.from_checkpoint, "refine-003")
        self.assertEqual(doctor.cycles, 3)
        listing = parser.parse_args(["checkpoints", "list", "run_001"])
        self.assertEqual(listing.checkpoint_action, "list")
        add = parser.parse_args([
            "checkpoints", "add", "run_001", "--name", "after coot",
            "--model", "fixed.pdb",
        ])
        self.assertEqual(add.name, "after coot")
        self.assertEqual(str(add.model), "fixed.pdb")
        use = parser.parse_args(["checkpoints", "use", "run_001", "refine-001"])
        self.assertEqual(use.checkpoint, "refine-001")
        active_use = parser.parse_args(["checkpoints", "use", "refine-001"])
        self.assertEqual(active_use.run_or_checkpoint, "refine-001")
        self.assertIsNone(active_use.checkpoint)
        show = parser.parse_args(["show", "last", "dataset", "--stage", "autosol"])
        self.assertEqual(show.target, "last")
        self.assertEqual(str(show.dataset), "dataset")
        self.assertEqual(show.stage, "autosol")
        inspect = parser.parse_args([
            "show", "run_001", "--checkpoint", "refine-005",
        ])
        self.assertEqual(inspect.checkpoint, "refine-005")


if __name__ == "__main__":
    unittest.main()
