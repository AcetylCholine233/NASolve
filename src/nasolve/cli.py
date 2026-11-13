"""NASolve command-line interface."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import os
import shlex
import sys
from pathlib import Path
from typing import Mapping

from .autorefine import AutoRefineError, execute_autorefine
from .autosol import AutoSolPreparationError, execute_autosol
from .autosol_metal import (
    DEFAULT_PARTIAL_MODEL_RMSD,
    FOUND as METAL_FOUND,
    execute_metal_search,
    parse_metals,
)
from .automr import AutoMRInputError, ModelAssessmentError, prepare_automr
from .checkpoints import (
    CheckpointError,
    add_checkpoint,
    list_checkpoints,
    select_checkpoint,
)
from .config import AppConfig, ConfigError, default_config_path, load_config, save_config
from .coot_runtime import (
    CootDiscoveryError,
    discover_coot,
    installation_from_candidate as coot_from_candidate,
    remember_coot,
)
from .coot_view import (
    CootViewError, launch_coot_files, launch_coot_view, resolve_run, resolve_view_profile,
)
from .backbone import BackboneError, requested_backbone_policy, write_backbone_review
from .phenix_runtime import (
    PhenixDiscoveryError,
    discover_phenix,
    installation_from_candidate,
    remember_phenix,
)
from .phaser import PhaserExecutionError, execute_phaser
from .postmr import PostMRPreparationError, prepare_postmr
from .refine_doctor import RefineDoctorError, execute_refine_doctor


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="nasolve", description="Guarded nucleic-acid structure solution")
    parser.add_argument("--phenix-root", help="one-run Phenix root, phenix_env.sh, or executable override")
    subparsers = parser.add_subparsers(dest="command", required=True)
    configure = subparsers.add_parser("configure", help="save external-tool configuration")
    configure_sub = configure.add_subparsers(dest="configure_target", required=True)
    phenix = configure_sub.add_parser("phenix", help="discover, validate, and remember Phenix")
    phenix.add_argument("path", nargs="?", help="Phenix root, phenix_env.sh, or phenix executable")
    coot = configure_sub.add_parser("coot", help="discover, validate, and remember Coot")
    coot.add_argument("path", nargs="?", help="Coot executable or installation directory")
    workspace = subparsers.add_parser(
        "workspace", help="remember or inspect the active dataset/run on this computer"
    )
    workspace_sub = workspace.add_subparsers(dest="workspace_action", required=True)
    workspace_use = workspace_sub.add_parser(
        "use", help="make a dataset or completed run the active target"
    )
    workspace_use.add_argument("target", type=Path)
    workspace_sub.add_parser("status", help="show the active dataset and run")
    workspace_sub.add_parser("clear", help="forget the active dataset and run")
    preset = subparsers.add_parser("preset", help="validate a versioned project preset")
    preset_sub = preset.add_subparsers(dest="preset_action", required=True)
    preset_check = preset_sub.add_parser("check", help="check a preset and its resource checksums")
    preset_check.add_argument("source", nargs="?", default="5w6w", help="5w6w or a preset TOML file")
    preset_check.add_argument("--json", action="store_true", help="print the structured preset record")
    campaign = subparsers.add_parser("campaign", help="plan, run, and inspect a frozen dataset campaign")
    campaign_sub = campaign.add_subparsers(dest="campaign_action", required=True)
    campaign_plan = campaign_sub.add_parser("plan", help="freeze input and model selection without running stages")
    campaign_plan.add_argument("root", nargs="?", type=Path, default=Path("."), help="parent directory containing datasets")
    campaign_plan.add_argument("--preset", default="5w6w", help="5w6w or a preset TOML file")
    campaign_plan.add_argument("--dataset", action="append", help="exact dataset subdirectory; repeat to select several")
    campaign_plan.add_argument("--frames-dir", type=Path, help="explicit approved MR_frames catalogue root")
    campaign_plan.add_argument("--json", action="store_true", help="print the complete structured plan")
    campaign_status = campaign_sub.add_parser("status", help="verify the existing plan and its frozen inputs")
    campaign_status.add_argument("root", nargs="?", type=Path, default=Path("."), help="campaign parent directory")
    campaign_status.add_argument("--json", action="store_true", help="print structured status and integrity issues")
    campaign_run = campaign_sub.add_parser("run", help="execute or resume selected datasets sequentially")
    campaign_run.add_argument("root", nargs="?", type=Path, default=Path("."), help="campaign parent directory")
    campaign_run.add_argument("--dataset", action="append", help="exact planned dataset name; repeat to select several")
    campaign_run.add_argument(
        "--through", choices=("preflight", "phaser", "postmr", "autosol", "autorefine"),
        default="autorefine", help="stop after this stage (default: autorefine)",
    )
    campaign_run.add_argument("--json", action="store_true", help="print structured final status without progress messages")
    campaign_pause = campaign_sub.add_parser("pause", help="request a pause after the active stage finishes")
    campaign_pause.add_argument("root", nargs="?", type=Path, default=Path("."), help="campaign parent directory")
    campaign_pause.add_argument("--json", action="store_true", help="print structured status")
    campaign_retry = campaign_sub.add_parser("retry", help="schedule a fresh attempt while preserving previous runs")
    campaign_retry.add_argument("root", nargs="?", type=Path, default=Path("."), help="campaign parent directory")
    campaign_retry.add_argument("--dataset", required=True, action="append", help="one exact planned dataset name")
    campaign_retry.add_argument("--json", action="store_true", help="print structured status")
    subparsers.add_parser("check", help="validate NASolve's runtime without solving a structure")
    automr = subparsers.add_parser(
        "automr", help="validate and freeze molecular-replacement inputs"
    )
    automr.add_argument(
        "dataset", nargs="?", type=Path,
        help="dataset directory (default: active workspace, then current directory)",
    )
    frames = automr.add_mutually_exclusive_group()
    frames.add_argument(
        "-W", "--w-frame", dest="frame", action="store_const", const="W",
        help="use the standard 5W6W/W frame",
    )
    frames.add_argument(
        "-3GBI", dest="frame", action="store_const", const="3GBI",
        help="use the standard 3GBI frame",
    )
    frames.add_argument("--frame", help="standard frame name (W, 5W6W, or 3GBI)")
    automr.add_argument("--pair", help="ordered standard-site pair, for example D:T")
    automr.add_argument(
        "--model",
        help=(
            "explicit MR model selector; in standard mode search the dataset and "
            "selected frame catalogue, otherwise use a dataset-relative model"
        ),
    )
    automr.add_argument(
        "--show", action="store_true",
        help="with --execute: open the Phaser model and map in Coot afterwards",
    )
    automr.add_argument(
        "--coot", help="with --show: one-run Coot executable override"
    )
    phaser_data = automr.add_mutually_exclusive_group()
    phaser_data.add_argument(
        "-I", dest="phaser_data", action="store_const", const="I",
        help="Phaser uses the MTZ's merged intensities (I/SIGI)",
    )
    phaser_data.add_argument(
        "-F", dest="phaser_data", action="store_const", const="F",
        help=(
            "Phaser uses the MTZ's merged amplitudes (F/SIGF); the default "
            "unless nasolve.txt sets phaser_data"
        ),
    )
    phaser_data.add_argument(
        "--phaser-data", dest="phaser_data", metavar="{F,I}",
        help="long form of -F/-I",
    )
    automr.add_argument(
        "--allow-p1-standard", action="store_true",
        help="strongly discouraged: allow a standard frame in P1 using three MR copies",
    )
    automr.add_argument(
        "--config", type=Path,
        help="input file (default: DATASET/nasolve.txt; generated there when absent)",
    )
    automr.add_argument(
        "--frames-dir", type=Path,
        help="approved MR frame directory (normally discovered from the repository/package)",
    )
    automr.add_argument(
        "--execute", action="store_true",
        help="run Phenix Phaser after the guarded preflight",
    )
    automr.add_argument(
        "--mirror", action="store_true",
        help="mirror the selected D/L nucleic-acid model with NARestraints before MR",
    )
    postmr = subparsers.add_parser(
        "postmr", help="prepare an accepted Phaser solution for refinement"
    )
    postmr.add_argument(
        "run", nargs="?", type=Path,
        help="completed AutoMR run directory (default: active workspace run)",
    )
    postmr.add_argument("--coot", help="one-run Coot executable override")
    postmr.add_argument(
        "--original-model", action="store_true",
        help=(
            "start from the original Phaser model even if metal sites were merged "
            "into it in Coot (Phaser/*-merged.pdb)"
        ),
    )
    postmr.add_argument(
        "--allow-mr-review", action="store_true",
        help="explicitly continue from a TFZ 7.0-7.99 MR_REVIEW result",
    )
    postmr.add_argument(
        "--modified-pairs-only",
        action="store_true",
        help="guess base pairs and restrain only pairs containing a modified nucleotide",
    )
    postmr.add_argument(
        "--allow-unreviewed-backbone", action="store_true",
        help=("experimental: continue at explicitly marked non-standard backbone "
              "sites without validating their linkage chemistry"),
    )
    backbone_review = subparsers.add_parser(
        "backbone-review",
        help="inspect experimental backbone sites in Coot and record human review",
    )
    backbone_review.add_argument(
        "run", nargs="?", type=Path,
        help="NASolve run (default: active workspace run)",
    )
    backbone_review.add_argument("--coot", help="one-run Coot executable override")
    autosol = subparsers.add_parser(
        "autosol",
        help="MR-SAD AutoSol: search for named metal ions after Phaser (--metal), "
             "or phase from a nucleotide heavy atom after PostMR",
    )
    autosol.add_argument(
        "run", nargs="?", type=Path,
        help="AutoMR run directory (default: active workspace run); without --metal it "
             "must be POSTMR_READY",
    )
    autosol.add_argument(
        "--metal", action="append", metavar="ELEMENT",
        help="search for this metal right after Phaser, e.g. --metal AG; repeat or use "
             "AG,HG for several (the first is AutoSol's atom_type)",
    )
    autosol.add_argument(
        "--rmsd", type=float, default=None, metavar="A",
        help="estimated RMSD of the Phaser model for the metal search "
             "(autosol.phasing.partpdb_rms; default 1.0)",
    )
    autosol.add_argument(
        "--sequence", type=Path, metavar="FILE",
        help="sequence file for the metal search (default: the frame's seq_base.txt "
             "frozen by AutoMR, else the chains of the Phaser model)",
    )
    autosol.add_argument(
        "--allow-mr-review", action="store_true",
        help="allow the metal search on a Phaser result that needs review (TFZ 7-8)",
    )
    autosol.add_argument(
        "--no-coot", action="store_true",
        help="do not open Coot when the metal search finds sites",
    )
    autosol.add_argument("--coot", help="one-run Coot executable override")
    autorefine = subparsers.add_parser(
        "autorefine", help="run one quiet, checkpointed Phenix refinement round"
    )
    autorefine.add_argument(
        "run", nargs="?", type=Path,
        help="completed PostMR/AutoSol run directory (default: active workspace run)",
    )
    autorefine.add_argument(
        "--from", dest="from_checkpoint",
        help="checkpoint ID or bookmark to branch from (default: current)",
    )
    autorefine.add_argument(
        "--recipe", default="AutoRefine/default",
        help="refinement recipe label recorded in provenance",
    )
    autorefine.add_argument(
        "--cycles", type=int, default=5,
        help="number of refinement macrocycles (default: 5)",
    )
    autofdp = subparsers.add_parser(
        "autofdp",
        help="measure metal occupancy by f'' refinement from a model and an MTZ "
             "(steps 1-2 of the plan implemented)",
        description=(
            "Hold coordinates, B and occupancy, fix f' at theory for the collection "
            "wavelength, refine f'' alone, and adjust occupancy until refined f'' is "
            "within the tolerance of theory. Needs only a starting model (PDB) and "
            "an MTZ with anomalous data and Free-R flags; labels, flags and "
            "wavelength are read from the MTZ, the metal sites from the model. Only "
            "steps 1 and 2 are implemented: copy the model and set each metal's B to "
            "the mean B of its neighbouring nucleotides (step 1), then run one "
            "phenix.refine round with every strategy off except the anomalous "
            "scatterers (f' fixed at theory, f'' refined from 0), open the result "
            "and its maps in Coot, and print each refined f'' against theory "
            "(step 2). Occupancy adjustment and the later steps are NOT IMPLEMENTED "
            "YET. See docs/autofdp.md."
        ),
    )
    autofdp.add_argument(
        "model", type=Path, metavar="MODEL.pdb",
        help="the starting model (PDB format), e.g. an AutoRefine refined_001.pdb",
    )
    autofdp.add_argument(
        "mtz", type=Path, metavar="DATA.mtz",
        help="the reflection file with anomalous data and Free-R flags "
             "(the two files may be given in either order)",
    )
    autofdp.add_argument(
        "--output", type=Path, metavar="DIR",
        help="folder in which a new run_NNN is made (default: AutoFDP beside the MTZ)",
    )
    autofdp.add_argument(
        "--site", action="append", default=[], metavar="SITE",
        help="metal site whose B step 1 sets, e.g. B:101; repeatable (default: every "
             "metal site). Every metal site is still refined in the f'' round",
    )
    autofdp.add_argument(
        "--wavelength", type=float, metavar="A",
        help="X-ray wavelength in A, when the MTZ header has none (default: MTZ header)",
    )
    autofdp.add_argument(
        "--cif", action="append", default=[], type=Path, metavar="FILE",
        help="restraint dictionary for a residue phenix.refine does not know; repeatable",
    )
    autofdp.add_argument(
        "--element", action="append", default=[], metavar="SITE=ELEMENT",
        help="(later step, not used yet) treat a site as another element, e.g. B:101=HG (renamed in a new checkpoint)",
    )
    autofdp.add_argument(
        "--candidates", metavar="ELEMENT[,ELEMENT...]",
        help="(later step, not used yet) elements compared in the identity check (default: AG,HG plus the model's)",
    )
    autofdp.add_argument(
        "--theory", action="append", default=[], metavar="ELEMENT=FDP",
        help="theoretical f'' to compare with, e.g. HG=10.17; repeatable "
             "(default: cctbx Henke at the MTZ wavelength)",
    )
    autofdp.add_argument(
        "--tolerance", type=float, default=1.0, metavar="E",
        help="refined f'' within E electrons of theory counts as matching (default: 1.0)",
    )
    autofdp.add_argument(
        "--max-iterations", type=int, default=8, metavar="N",
        help="(later step, not used yet) f'' rounds before review (default: 8)",
    )
    autofdp.add_argument(
        "--cycles", type=int, default=3, metavar="N",
        help="phenix.refine macrocycles in the f'' round (default: 3)",
    )
    autofdp.add_argument(
        "--b-start", default="surrounding", metavar="surrounding|model|VALUE",
        help="site B at the start: mean of coordinating nucleotides (default), the model's, or a number",
    )
    autofdp.add_argument(
        "--b-step", type=float, default=20.0, metavar="A2",
        help="(later step, not used yet) B change when the difference map calls for B instead of occupancy (default: 20)",
    )
    autofdp.add_argument(
        "--occupancy-floor", type=float, default=0.05, metavar="Q",
        help="(later step, not used yet) below this occupancy the site goes to identity review (default: 0.05)",
    )
    autofdp.add_argument(
        "--map-sigma", type=float, default=3.0, metavar="S",
        help="(later step, not used yet) mFo-DFc peak or hole at the site that counts as significant (default: 3.0)",
    )
    autofdp.add_argument(
        "--no-confirm-b", action="store_true",
        help="(later step, not used yet) skip the round that refines only the sites' B at fixed occupancy and f''",
    )
    autofdp.add_argument(
        "--final", action="store_true",
        help="(later step, not used yet) after convergence, refine XYZ, B, occupancy and f'' together (f' stays fixed)",
    )
    autofdp.add_argument(
        "--plan", action="store_true",
        help="check the inputs and show what was read (labels, Free-R flags, "
             "wavelength, theory, B changes) without refining or writing anything",
    )
    autofdp.add_argument(
        "--no-coot", action="store_true",
        help="do not open the f'' round's model and maps in Coot",
    )
    autofdp.add_argument(
        "--coot", type=Path, metavar="PATH",
        help="Coot executable or installation (default: the configured Coot)",
    )
    autofdp.add_argument(
        "--phenix-root", type=Path, metavar="PATH",
        help="Phenix installation (default: the configured Phenix)",
    )
    doctor = subparsers.add_parser(
        "refine-doctor",
        aliases=["refine_doctor"],
        help="audit and triage a refinement with bounded checkpoint branches",
    )
    doctor.add_argument(
        "run", nargs="?", type=Path,
        help="completed PostMR/AutoRefine run directory (default: active workspace run)",
    )
    doctor.add_argument(
        "--from", dest="from_checkpoint",
        help="checkpoint ID or bookmark to diagnose (default: current)",
    )
    doctor.add_argument(
        "--cycles", type=int, default=3,
        help="macrocycles per bounded diagnostic branch (default: 3)",
    )
    doctor.add_argument(
        "--max-trials", type=int, default=5,
        help="maximum eligible Doctor trials, 1-10 (default: 5)",
    )
    checkpoints = subparsers.add_parser(
        "checkpoints", help="list, bookmark, import, or select refinement checkpoints"
    )
    checkpoint_sub = checkpoints.add_subparsers(dest="checkpoint_action", required=True)
    checkpoint_list = checkpoint_sub.add_parser("list", help="list the checkpoint tree")
    checkpoint_list.add_argument("run", nargs="?", type=Path)
    checkpoint_add = checkpoint_sub.add_parser(
        "add", help="bookmark the current checkpoint or import a manual model"
    )
    checkpoint_add.add_argument("run", nargs="?", type=Path)
    checkpoint_add.add_argument("--name", required=True, help="user-visible checkpoint name")
    checkpoint_add.add_argument("--model", type=Path, help="manual PDB model to import")
    checkpoint_add.add_argument(
        "--mtz", type=Path,
        help="deliberate replacement observation MTZ for an imported model",
    )
    checkpoint_add.add_argument(
        "--from", dest="from_checkpoint",
        help="parent checkpoint ID or bookmark (default: current)",
    )
    checkpoint_use = checkpoint_sub.add_parser("use", help="select a reusable checkpoint")
    checkpoint_use.add_argument(
        "run_or_checkpoint",
        help="run path, or checkpoint ID/bookmark when an active run is selected",
    )
    checkpoint_use.add_argument(
        "checkpoint", nargs="?",
        help="checkpoint ID/bookmark when an explicit run path is supplied",
    )
    show = subparsers.add_parser(
        "show", help="open a completed run stage in an isolated graphical Coot session"
    )
    show.add_argument(
        "target", nargs="?",
        help=(
            "AutoMR/run_NNN or AutoSolMetal/run_NNN folder (relative paths are "
            "tried from here, then from the active dataset), 'last', or omit to "
            "use the active workspace run"
        ),
    )
    show.add_argument(
        "dataset", nargs="?", type=Path,
        help="dataset directory required by 'show last'",
    )
    show.add_argument(
        "--stage", choices=("automr", "metal", "postmr", "autosol", "autorefine"),
        help="view a specific stage instead of the selected current checkpoint",
    )
    show.add_argument(
        "--checkpoint",
        help="open a specific checkpoint without changing the current selection",
    )
    show.add_argument("--coot", help="one-run Coot executable override")
    return parser


def _narestraints_version() -> str:
    try:
        return importlib.metadata.version("NARestraints")
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def _run_report_path(path: Path) -> Path:
    selected = path.expanduser().resolve()
    return selected if selected.name == "report.json" else selected / "report.json"


def _workspace_dataset(
    explicit: Path | None,
    config: AppConfig,
    *,
    current_directory_fallback: bool = True,
) -> Path:
    value: Path | None = explicit
    if value is None and config.workspace.dataset:
        value = Path(config.workspace.dataset)
    if value is None and current_directory_fallback:
        value = Path(".")
    if value is None:
        raise ConfigError(
            "No active dataset. Run 'nasolve workspace use DATASET' or pass a path."
        )
    dataset = value.expanduser().resolve()
    if not dataset.is_dir():
        raise ConfigError(f"Workspace dataset does not exist: {dataset}")
    return dataset


def _workspace_run(explicit: Path | None, config: AppConfig) -> Path:
    value = explicit or (Path(config.workspace.run) if config.workspace.run else None)
    if value is None:
        raise ConfigError(
            "No active run. Run 'nasolve workspace use RUN' or pass a run path."
        )
    report_path = _run_report_path(value)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigError(f"Active run has no readable report: {report_path}") from exc
    if not isinstance(report, Mapping):
        raise ConfigError(f"Active run report is not a JSON object: {report_path}")
    if report.get("workflow") != "automr":
        raise ConfigError(f"Active run is not a NASolve AutoMR run: {report_path.parent}")
    return report_path.parent


def _remember_workspace(config: AppConfig, run: Path) -> None:
    resolved = run.expanduser().resolve()
    config.workspace.run = str(resolved)
    config.workspace.dataset = str(
        resolved.parent.parent if resolved.parent.name == "AutoMR" else resolved.parent
    )


def _workspace(args: argparse.Namespace) -> int:
    try:
        config = load_config()
        if args.workspace_action == "clear":
            config.workspace.dataset = None
            config.workspace.run = None
            saved_at = save_config(config)
            print("Active NASolve workspace cleared.")
            print(f"Saved: {saved_at}")
            return 0
        if args.workspace_action == "use":
            target = args.target.expanduser().resolve()
            report_path = _run_report_path(target)
            if report_path.is_file():
                run = _workspace_run(target, config)
                _remember_workspace(config, run)
            else:
                if not target.is_dir():
                    raise ConfigError(f"Workspace target does not exist: {target}")
                config.workspace.dataset = str(target)
                config.workspace.run = None
            saved_at = save_config(config)
            print(f"Active dataset: {config.workspace.dataset}")
            print(f"Active run: {config.workspace.run or 'none'}")
            print(f"Saved: {saved_at}")
            return 0
        print(f"Active dataset: {config.workspace.dataset or 'none'}")
        print(f"Active run: {config.workspace.run or 'none'}")
        print(f"Configuration: {default_config_path()}")
        return 0
    except ConfigError as exc:
        print(f"Workspace error: {exc}", file=sys.stderr)
        return 2


def _preset(args: argparse.Namespace) -> int:
    from .presets import PresetError, load_preset

    try:
        preset = load_preset(args.source)
    except PresetError as exc:
        print(f"Preset error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(preset.to_dict(), indent=2, sort_keys=True))
    else:
        print(f"Preset: {preset.id} {preset.version}")
        print(f"Source: {preset.source}")
        print(f"Configuration SHA-256: {preset.config_sha256}")
        print(f"Validated resources: {len(preset.resources)}")
        print("Recipe 5'-phosphate/OP3 sites: " + (", ".join(preset.terminal_phosphate_sites) or "none"))
    return 0


def _campaign(args: argparse.Namespace) -> int:
    from .campaigns import CampaignError, plan_campaign

    try:
        if args.campaign_action == "plan":
            result = plan_campaign(
                args.root,
                preset=args.preset,
                datasets=tuple(args.dataset) if args.dataset is not None else None,
                frames_directory=args.frames_dir,
            )
        else:
            from .campaign_execution import (
                execute_campaign, execution_status, request_pause, retry_dataset,
            )

            if args.campaign_action == "run":
                result = execute_campaign(
                    args.root,
                    datasets=tuple(args.dataset) if args.dataset is not None else None,
                    through=args.through,
                    phenix_root=args.phenix_root,
                    progress=None if args.json else lambda message: print(message, flush=True),
                )
            elif args.campaign_action == "pause":
                result = request_pause(args.root)
            elif args.campaign_action == "retry":
                if len(args.dataset) != 1:
                    raise CampaignError("campaign retry requires exactly one --dataset")
                result = retry_dataset(args.root, args.dataset[0])
            else:
                result = execution_status(args.root)
    except CampaignError as exc:
        print(f"Campaign error: {exc}", file=sys.stderr)
        return 2
    counts = result["counts"]
    execution = result.get("execution")
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"Campaign: {args.root.expanduser().resolve()}")
        print(f"Plan: {result['plan_path']}")
        print(f"Preset: {result['preset']['id']} {result['preset']['version']}")
        if execution is None:
            print("Input and model selection only; scientific preflight and stage execution are pending.")
            print(
                f"Datasets: {counts['total']} total, {counts['discovered']} discovered, "
                f"{counts['blocked']} blocked"
            )
        else:
            print(f"Execution: {execution['state']}")
            print(f"Datasets: {counts['total']} planned")
            if execution.get("pause_requested"):
                print("Pause requested; the active stage will finish before execution stops.")
        if "integrity" in result:
            print(f"Frozen input integrity: {result['integrity']}")
        progress_by_id = {
            item["id"]: item for item in execution["datasets"]
        } if execution is not None else {}
        for planned in result["datasets"]:
            item = {**planned, **progress_by_id.get(planned["id"], {})}
            integrity = f"; integrity {item['integrity']}" if "integrity" in item else ""
            print(f"  {item['id']}: {item['status']}{integrity}")
            effective = planned.get("effective_config")
            if effective is not None:
                from .phosphate import phosphate_intent_summary
                print("    " + phosphate_intent_summary(
                    effective.get("allow_op3_sites", []), effective.get("phosphate_intent")))
                if effective.get("sequence_reference"):
                    print(f"    Sequence reference: {effective['sequence_reference']}")
                if effective.get("sequence_thread"):
                    print(f"    Sequence thread: {effective['sequence_thread']['id']}")
                if effective.get("model_selector"):
                    provider = effective.get("model_provider") or {}
                    location = provider.get("location", "explicit")
                    print(
                        f"    Model override: {effective['model_selector']} "
                        f"({location})"
                    )
            if item.get("diagnostic"):
                print(f"    {item['diagnostic']}")
            if item.get("run"):
                print(f"    Run: {args.root.expanduser().resolve() / item['run']}")
            if item.get("next_stage"):
                print(f"    Next stage: {item['next_stage']}")
            if item.get("checkpoint"):
                print(f"    Checkpoint: {item['checkpoint']}")
                if item.get("run"):
                    command = [
                        "./nasolve", "show", str(args.root.expanduser().resolve() / item["run"]),
                        "--checkpoint", item["checkpoint"],
                    ]
                    print(f"    Inspect: {shlex.join(command)}")
            if item.get("numerical_success"):
                print("    Numerical refinement criteria passed; inspect the model and maps.")
            if item.get("duplicate_of"):
                print(f"    Same observation bytes as {item['duplicate_of']}")
            for issue in item.get("integrity_issues", []):
                print(f"    {issue}")
        for issue in result.get("integrity_issues", []):
            print(f"  {issue}")
    flagged = counts["blocked"] or result.get("integrity") == "DRIFT"
    if execution is not None:
        flagged = flagged or any(
            item["status"] in {"BLOCKED", "NO_SOLUTION", "AWAITING_INSPECTION", "DRIFT"}
            for item in execution["datasets"]
        )
    return 3 if flagged else 0


def _configure_phenix(path: str | None) -> int:
    candidate = path
    if not candidate:
        if not sys.stdin.isatty():
            print("Phenix path is required in non-interactive mode.", file=sys.stderr)
            return 2
        candidate = input("Phenix installation folder or phenix_env.sh: ").strip()
    try:
        installation = installation_from_candidate(candidate)
        config = load_config()
        remember_phenix(config, installation)
        saved_at = save_config(config)
    except (ConfigError, PhenixDiscoveryError) as exc:
        print(f"Phenix configuration failed: {exc}", file=sys.stderr)
        return 1
    print(f"Phenix {installation.version}: OK")
    print(f"Source: {installation.source}")
    print(f"Root: {installation.root}")
    print(f"Saved: {saved_at}")
    return 0


def _configure_coot(path: str | None) -> int:
    candidate = path
    if not candidate:
        if not sys.stdin.isatty():
            print("Coot path is required in non-interactive mode.", file=sys.stderr)
            return 2
        candidate = input("Coot executable or installation folder: ").strip()
    try:
        installation = coot_from_candidate(candidate)
        config = load_config()
        remember_coot(config, installation)
        saved_at = save_config(config)
    except (ConfigError, CootDiscoveryError) as exc:
        print(f"Coot configuration failed: {exc}", file=sys.stderr)
        return 1
    print(f"Coot {installation.version}: OK")
    print(f"Executable: {installation.executable}")
    print(f"Saved: {saved_at}")
    return 0


def _check(explicit: str | None) -> int:
    print(f"Python {sys.version.split()[0]}: OK")
    na_version = _narestraints_version()
    print(f"NARestraints {na_version}: " + ("OK" if na_version != "not installed" else "MISSING"))
    try:
        config = load_config()
        installation = discover_phenix(config, explicit=explicit)
        saved_at = None
        if explicit is None:
            remember_phenix(config, installation)
            saved_at = save_config(config)
    except (ConfigError, PhenixDiscoveryError) as exc:
        print(f"Phenix: NOT CONFIGURED\n{exc}")
        print(f"Configuration file: {default_config_path()}")
        return 1
    print(f"Phenix {installation.version}: OK")
    print(f"Source: {installation.source}")
    for name, executable in installation.executables.items():
        print(f"  {name}: {executable}")
    if saved_at:
        print(f"Remembered: {saved_at}")
    try:
        coot = discover_coot(config)
        remember_coot(config, coot)
        save_config(config)
    except CootDiscoveryError as exc:
        print("Coot: NOT CONFIGURED (needed only when PostMR must mutate supported bases)")
        print(str(exc))
    else:
        print(f"Coot {coot.version}: OK")
        print(f"  coot: {coot.executable}")
    return 0 if na_version != "not installed" else 1


def _print_metal_merge_hint(merge) -> None:
    if merge is None:
        return
    print(
        "To keep the sites: in Coot click the NASolve toolbar button -> "
        "'Merge metal sites and save' (or run nasolve_merge_metal_sites() in "
        "Calculate > Scripting > Python)."
    )
    print(f"Merged model will be saved as: {merge.output_directory / (merge.merged_stem + '.pdb')}")
    print(f"'nasolve postmr {merge.run_directory}' then starts from the merged model.")
    for warning in merge.warnings:
        print(f"Warning: {warning}")


def _open_phaser_view(args: argparse.Namespace, phaser) -> None:
    """Open Phaser's model and map in Coot; a Coot problem never fails the run."""
    later = f"Open in Coot with: nasolve show {phaser.run_directory} --stage automr"
    if phaser.solution_pdb is None or phaser.solution_mtz is None:
        return
    if not args.show:
        print(later)
        return
    try:
        config = load_config()
        coot = discover_coot(config, explicit=args.coot)
        remember_coot(config, coot)
        save_config(config)
        view = launch_coot_view(phaser.run_directory, coot.executable, stage="automr")
    except (ConfigError, CootDiscoveryError, CootViewError) as exc:
        print(f"Coot view error: {exc}", file=sys.stderr)
        print(later)
        return
    print(f"Opened Coot (PID {view.pid}) with the Phaser model and map")
    print(f"Coot log: {view.log_path}")


def _automr(args: argparse.Namespace) -> int:
    if args.show and not args.execute:
        print("AutoMR input error: --show applies only with --execute", file=sys.stderr)
        return 2
    if args.coot and not args.show:
        print("AutoMR input error: --coot applies only with --show", file=sys.stderr)
        return 2
    try:
        config = load_config()
        dataset = _workspace_dataset(args.dataset, config)
        installation = discover_phenix(config, explicit=args.phenix_root)
        if args.phenix_root is None:
            remember_phenix(config, installation)
            save_config(config)
        result = prepare_automr(
            dataset,
            config_path=args.config,
            frame_override=args.frame,
            pair_override=args.pair,
            model_override=args.model,
            frames_dir=args.frames_dir,
            allow_p1_standard=args.allow_p1_standard,
            mirror=args.mirror,
            mtz_dump_executable=installation.executables.get("phenix.mtz.dump"),
            phenix_environment=installation.environment,
            phaser_data=args.phaser_data,
        )
        _remember_workspace(config, result.run_directory)
        save_config(config)
        if args.execute:
            phaser = execute_phaser(
                result.report_path,
                installation.executables["phenix.phaser"],
                environment=installation.environment,
                phenix_version=installation.version,
            )
    except (
        AutoMRInputError,
        ModelAssessmentError,
        ConfigError,
        PhenixDiscoveryError,
        PhaserExecutionError,
    ) as exc:
        print(f"AutoMR input error: {exc}", file=sys.stderr)
        return 2
    if args.execute:
        print(f"{phaser.status}: {phaser.message}")
        print("Phaser executed: yes")
        print(f"Run directory: {phaser.run_directory}")
        print(f"Phaser log: {phaser.log_path}")
        print(f"Report: {phaser.report_path}")
        if result.phosphate_summary:
            print(result.phosphate_summary)
        if phaser.solution_pdb:
            print(f"MR model: {phaser.solution_pdb}")
        if phaser.solution_mtz:
            print(f"MR reflections: {phaser.solution_mtz}")
        _open_phaser_view(args, phaser)
        return phaser.exit_code
    print(f"{result.status}: {result.message}")
    print("Phaser executed: no (preflight milestone)")
    print(f"Run directory: {result.run_directory}")
    print(f"Report: {result.report_path}")
    if result.phosphate_summary:
        print(result.phosphate_summary)
    if result.generated_config:
        print(f"Generated: {dataset / 'nasolve.txt'}")
    return 0


def _postmr(args: argparse.Namespace) -> int:
    try:
        config = load_config()
        run = _workspace_run(args.run, config)
        phenix = discover_phenix(config, explicit=args.phenix_root)
        if args.phenix_root is None:
            remember_phenix(config, phenix)
        coot = None
        try:
            coot = discover_coot(config, explicit=args.coot)
        except CootDiscoveryError:
            if args.coot:
                raise
        else:
            remember_coot(config, coot)
        save_config(config)
        result = prepare_postmr(
            run,
            phenix.executables["phenix.ready_set"],
            coot_executable=coot.executable if coot else None,
            environment=phenix.environment,
            allow_mr_review=args.allow_mr_review,
            modified_pairs_only=args.modified_pairs_only,
            allow_unreviewed_backbone=args.allow_unreviewed_backbone,
            use_merged_model=not args.original_model,
        )
    except (
        ConfigError,
        CootDiscoveryError,
        PhenixDiscoveryError,
        PostMRPreparationError,
    ) as exc:
        print(f"PostMR preparation error: {exc}", file=sys.stderr)
        return 2
    print(f"{result.status}: {result.message}")
    print(f"Run directory: {result.run_directory}")
    if result.input_model is not None:
        source = (
            "Phaser model with merged metal sites"
            if result.input_model_source == "metal-merged" else "Phaser model"
        )
        print(f"Starting model: {result.input_model} ({source})")
    print(f"Prepared model: {result.model_path}")
    for site in result.metal_sites:
        print(
            f"Anomalous scatterer: {site['element']} {site['site']} {site['atom_name']} "
            "(AutoRefine refines its f' and f'')"
        )
    if result.anomalous_scatterer_phil is not None:
        print(f"Anomalous-scatterer PHIL for phenix.refine by hand: {result.anomalous_scatterer_phil}")
    print(f"ReadySet log: {result.readyset_log}")
    for restraint in result.restraint_paths:
        print(f"Restraint: {restraint}")
    print(f"Report: {result.report_path}")
    try:
        run_report = json.loads((result.run_directory / "report.json").read_text(encoding="utf-8"))
        policy = requested_backbone_policy(run_report)
        pending = tuple(policy["experimental_passthrough_sites"])
    except (OSError, ValueError, BackboneError):
        pending = ()
    if pending:
        print(_color("WARNING: unreviewed non-standard backbone chemistry: " + ", ".join(pending), "33"))
        print("NASolve did not validate those linkage atoms. After refinement, run:")
        print("  " + shlex.join(["./nasolve", "backbone-review", str(result.run_directory)]))
    return 0


def _backbone_review(args: argparse.Namespace) -> int:
    try:
        config = load_config()
        run = _workspace_run(args.run, config)
        report = json.loads((run / "report.json").read_text(encoding="utf-8"))
        policy = requested_backbone_policy(report)
        sites = tuple(policy["experimental_passthrough_sites"])
        if not sites:
            print("No experimental backbone passthrough sites are recorded for this run.")
            return 0
        if not sys.stdin.isatty():
            raise ConfigError("backbone-review is interactive; run it in a terminal")
        print(_color("NON-STANDARD BACKBONE CHEMISTRY — REVIEW REQUIRED", "33"))
        print("Flagged site(s): " + ", ".join(sites))
        print("NASolve intentionally did not infer or validate their inter-residue linkage chemistry.")
        answer = input("Show the flagged result in Coot? [y/N]: " ).strip().casefold()
        if answer not in {"y", "yes"}:
            print("Review remains pending.")
            return 2
        coot = discover_coot(config, explicit=args.coot)
        remember_coot(config, coot)
        save_config(config)
        view = launch_coot_view(run, coot.executable)
        print(f"Opened {view.model_path} in Coot (PID {view.pid}).")
        print("Inspect the flagged backbone site(s) and maps, then return here.")
        confirm = input("Confirm that you reviewed the non-standard backbone chemistry? [y/N]: " ).strip().casefold()
        if confirm not in {"y", "yes"}:
            print("Review remains pending; passthrough provenance is unchanged.")
            return 2
        review_dir = run / "PostMR" / "BackboneReview"
        review_dir.mkdir(parents=True, exist_ok=True)
        number = 1
        while (review_dir / f"review_{number:03d}.json").exists():
            number += 1
        profile = resolve_view_profile(run)
        review_path = review_dir / f"review_{number:03d}.json"
        write_backbone_review(
            review_path, sites=sites, model=profile.model_path, run=run, reviewed=True
        )
        print(_color("Backbone review recorded as USER_REVIEWED.", "32"))
        print(f"Review record: {review_path}")
        print("Experimental passthrough provenance remains in the run report.")
        return 0
    except (ConfigError, CootDiscoveryError, CootViewError, BackboneError, OSError, ValueError) as exc:
        print(f"Backbone review error: {exc}", file=sys.stderr)
        return 2


def _autosol_metal(args: argparse.Namespace) -> int:
    """Metal-ion search right after Phaser; opens Coot when sites are found."""
    try:
        metals = parse_metals(args.metal)
        config = load_config()
        run = _workspace_run(args.run, config)
        phenix = discover_phenix(config, explicit=args.phenix_root)
        if args.phenix_root is None:
            remember_phenix(config, phenix)
            save_config(config)
        autosol_executable = phenix.executables.get("phenix.autosol")
        mtz_dump_executable = phenix.executables.get("phenix.mtz.dump")
        if autosol_executable is None:
            raise AutoSolPreparationError(
                "The discovered Phenix installation has no phenix.autosol executable"
            )
        if mtz_dump_executable is None:
            raise AutoSolPreparationError("The metal search requires phenix.mtz.dump")
        result = execute_metal_search(
            run,
            autosol_executable,
            mtz_dump_executable,
            metals,
            partial_model_rmsd=(
                DEFAULT_PARTIAL_MODEL_RMSD if args.rmsd is None else args.rmsd
            ),
            allow_mr_review=args.allow_mr_review,
            sequence_file=args.sequence,
            environment=phenix.environment,
        )
    except (AutoSolPreparationError, ConfigError, PhenixDiscoveryError) as exc:
        print(f"AutoSol metal search error: {exc}", file=sys.stderr)
        return 2
    print(f"{result.status}: {result.message}")
    print(f"Metals searched: {', '.join(result.metals)}")
    print(f"Partial-model RMSD: {result.partial_model_rmsd:g} A")
    try:
        record = json.loads(result.report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        record = {}
    if record.get("sequence_source"):
        print(f"Sequence: {record['sequence_source']}")
    for warning in record.get("sequence_warnings") or ():
        print(f"Sequence warning: {warning}")
    print(f"AutoMR run: {result.run_directory}")
    print(f"Metal-search folder: {result.output_directory}")
    print(f"AutoSol log: {result.log_path}")
    for site in result.sites:
        near = ", ".join(
            f"{item['chain']}:{item['residue_number']} {item['residue_name']} "
            f"{item['atom_name']} {float(item['distance']):.2f} A"
            for item in site["nearest_residues"]  # type: ignore[union-attr]
        )
        print(
            f"  site {site['index']}: {site['element']} occupancy "
            f"{float(site['occupancy']):.2f}; nearest {near or 'no model atoms'}"
        )
    print(f"Report: {result.report_path}")
    if result.status != METAL_FOUND:
        print("Continue with 'nasolve postmr' as usual; nothing was added to the model.")
        return 0
    if result.heavy_atom_model is not None:
        print(f"Heavy-atom model: {result.heavy_atom_model}")
    print(f"Phaser model: {result.phaser_model}")
    print(f"Phaser map: {result.phaser_map}")
    later = f"Open later with: nasolve show {result.output_directory}"
    if args.no_coot:
        print(later)
        return 0
    try:
        coot = discover_coot(config, explicit=args.coot)
        remember_coot(config, coot)
        save_config(config)
        view = launch_coot_view(result.output_directory, coot.executable, stage="metal")
    except (ConfigError, CootDiscoveryError, CootViewError) as exc:
        print(f"Coot view error: {exc}", file=sys.stderr)
        print(later)
        return 0
    print(f"Opened Coot (PID {view.pid}) with the Phaser model, AutoSol sites and Phaser map")
    print(f"Coot log: {view.log_path}")
    _print_metal_merge_hint(view.metal_merge)
    return 0


def _autosol(args: argparse.Namespace) -> int:
    if args.metal:
        return _autosol_metal(args)
    if args.rmsd is not None or args.allow_mr_review or args.no_coot or args.sequence:
        print(
            "AutoSol error: --rmsd, --allow-mr-review, --sequence and --no-coot apply "
            "only to the metal search; add --metal ELEMENT",
            file=sys.stderr,
        )
        return 2
    try:
        config = load_config()
        run = _workspace_run(args.run, config)
        phenix = discover_phenix(config, explicit=args.phenix_root)
        if args.phenix_root is None:
            remember_phenix(config, phenix)
            save_config(config)
        autosol_executable = phenix.executables.get("phenix.autosol")
        mtz_dump_executable = phenix.executables.get("phenix.mtz.dump")
        if autosol_executable is None:
            raise AutoSolPreparationError(
                "The discovered Phenix installation has no phenix.autosol executable"
            )
        if mtz_dump_executable is None:
            raise AutoSolPreparationError(
                "AutoSol validation requires phenix.mtz.dump"
            )
        result = execute_autosol(
            run,
            autosol_executable,
            mtz_dump_executable,
            environment=phenix.environment,
        )
    except (
        AutoSolPreparationError,
        ConfigError,
        PhenixDiscoveryError,
    ) as exc:
        print(f"AutoSol preparation error: {exc}", file=sys.stderr)
        return 2
    print(f"{result.status}: {result.message}")
    print(f"Run directory: {result.run_directory}")
    print(f"AutoSol log: {result.log_path}")
    if result.heavy_atom_model:
        print(f"Heavy-atom model: {result.heavy_atom_model}")
    if result.refinement_data:
        print(f"Refinement data: {result.refinement_data}")
    if result.matched_distance is not None:
        print(f"Matched HA distance: {result.matched_distance:.3f} A")
    if result.status != "AUTOSOL_READY":
        print("Refinement continuation: yes, without automatically approved AutoSol phases")
    print(f"Report: {result.report_path}")
    return 0


def _color(text: str, code: str) -> str:
    if not sys.stdout.isatty() or "NO_COLOR" in os.environ:
        return text
    return f"\033[{code}m{text}\033[0m"


def _status_text(status: str, *, current: bool = False) -> str:
    colors = {
        "READY": "32",
        "SUCCESS": "32",
        "USER_APPROVED": "32",
        "REVIEW": "33",
        "FAILED": "31",
    }
    marker = _color(">", "36") if current else " "
    return f"{marker} {_color(status, colors.get(status, '0'))}"


def _format_metric(value: object, digits: int = 3) -> str:
    return f"{value:.{digits}f}" if isinstance(value, (int, float)) else "—"


def _autorefine(args: argparse.Namespace) -> int:
    try:
        config = load_config()
        run = _workspace_run(args.run, config)
        phenix = discover_phenix(config, explicit=args.phenix_root)
        if args.phenix_root is None:
            remember_phenix(config, phenix)
            save_config(config)
        mtz_dump = phenix.executables.get("phenix.mtz.dump")
        if mtz_dump is None:
            raise AutoRefineError("AutoRefine requires phenix.mtz.dump")

        def progress(checkpoint_id: str, log_path: Path) -> None:
            print(f"AutoRefine {checkpoint_id} started; full output: {log_path}")

        result = execute_autorefine(
            run,
            phenix.executables["phenix.refine"],
            mtz_dump,
            phenix_version=phenix.version,
            environment=phenix.environment,
            from_checkpoint=args.from_checkpoint,
            recipe=args.recipe,
            macro_cycles=args.cycles,
            progress=progress,
        )
    except (AutoRefineError, ConfigError, PhenixDiscoveryError) as exc:
        print(f"AutoRefine error: {exc}", file=sys.stderr)
        return 2
    status_color = (
        "32" if result.status == "AUTOREFINE_READY"
        else "33" if result.status in {"AUTOREFINE_REVIEW", "AUTOREFINE_ANOMALOUS_FALLBACK"}
        else "31"
    )
    print(f"\n{_color(result.status, status_color)}: {result.message}")
    print(f"Checkpoint: {result.checkpoint_id} (parent: {result.parent_checkpoint})")
    print("Cycles: " + str(args.cycles))
    payload = json.loads(result.report_path.read_text(encoding="utf-8"))
    inputs = payload.get("inputs", {})
    refinement = payload.get("refinement", {})
    if isinstance(inputs, dict) and isinstance(refinement, dict):
        print("Target: Automatic")
        print("Observations: " + ",".join(inputs.get("observation_labels", [])))
        phase_labels = inputs.get("phase_labels", [])
        print("Phases: " + (",".join(phase_labels) if phase_labels else "none"))
        selections = refinement.get("anomalous_selections", [])
        print("Anomalous: " + ("; ".join(selections) if selections else "off"))
    stats = result.statistics
    print("\n                 Initial    Final")
    print(
        "Rwork           "
        f"{_format_metric(stats.get('initial_r_work')):>9}"
        f"{_format_metric(stats.get('r_work')):>9}"
    )
    print(
        "Rfree           "
        f"{_format_metric(stats.get('initial_r_free')):>9}"
        f"{_format_metric(stats.get('r_free')):>9}"
    )
    print(f"Rfree - Rwork   {'':>9}{_format_metric(stats.get('r_free_minus_r_work')):>9}")
    print(f"Clashscore      {'':>9}{_format_metric(stats.get('clashscore'), 2):>9}")
    bond = _format_metric(stats.get("bond_rmsd"), 3)
    angle = _format_metric(stats.get("angle_rmsd"), 2)
    print(f"Bond RMSD       {'':>9}{(bond + ' A') if bond != '—' else bond:>9}")
    print(f"Angle RMSD      {'':>9}{(angle + ' deg') if angle != '—' else angle:>9}")
    anomalous_stats = stats.get("anomalous_scatterers")
    if isinstance(anomalous_stats, list) and anomalous_stats:
        print("\nAnomalous-scatterer strength:")
        for item in anomalous_stats:
            if not isinstance(item, dict):
                continue
            label = f"{item.get('site')} {item.get('atom_name')} ({item.get('element')})"
            wavelength = _format_metric(item.get("wavelength"), 6)
            refined_fdp = _format_metric(item.get("refined_f_double_prime"), 3)
            calculated_fdp = _format_metric(item.get("calculated_f_double_prime"), 3)
            occupancy = _format_metric(item.get("model_occupancy"), 3)
            apparent = _format_metric(item.get("apparent_anomalous_occupancy"), 3)
            b_factor = _format_metric(item.get("b_factor"), 2)
            expected_at_occ = _format_metric(
                item.get("calculated_f_double_prime_at_occupancy"), 3
            )
            source = item.get("wavelength_source")
            origin = f" (from {source})" if isinstance(source, str) and wavelength != "—" else ""
            print(f"  {label}: wavelength {wavelength} A{origin}")
            print(f"    f'' refined {refined_fdp}; calculated {calculated_fdp} e-")
            print(
                f"    occupancy {occupancy}; calculated f'' x occupancy "
                f"{expected_at_occ}; apparent anomalous occupancy {apparent}"
            )
            resolution = item.get("resolution_limit")
            if isinstance(resolution, (int, float)):
                expected_edge = _format_metric(
                    item.get("calculated_contribution_at_resolution_limit"), 3
                )
                refined_edge = _format_metric(
                    item.get("refined_contribution_at_resolution_limit"), 3
                )
                print(
                    f"    B {b_factor} A^2; at {resolution:.3f} A: "
                    f"calculated {expected_edge}, refined {refined_edge} e-"
                )
    series = stats.get("cycle_series")
    if isinstance(series, list) and series:
        work_trend = " -> ".join(
            _format_metric(item.get("r_work"))
            for item in series if isinstance(item, dict)
        )
        free_trend = " -> ".join(
            _format_metric(item.get("r_free"))
            for item in series if isinstance(item, dict)
        )
        print(f"\nRwork trend: {work_trend}")
        print(f"Rfree trend: {free_trend}")
    diagnostics = stats.get("diagnostics")
    if isinstance(diagnostics, list) and diagnostics:
        print("\nActionable Phenix messages:")
        for message in diagnostics[:5]:
            print(f"  - {message}")
        if len(diagnostics) > 5:
            print(f"  - {len(diagnostics) - 5} more in the full log")
    if result.model_path:
        print(f"\nRefined model: {result.model_path}")
    if result.model_cif:
        print(f"Refined model mmCIF: {result.model_cif}")
    if result.reflection_cif:
        print(f"Reflections mmCIF: {result.reflection_cif}")
    if result.map_coefficients:
        print(f"Map coefficients: {result.map_coefficients}")
    print(f"Full log: {result.log_path}")
    print(f"Report: {result.report_path}")
    if result.selected_as_current:
        print(_color("Current checkpoint updated.", "32"))
    else:
        print("Current checkpoint unchanged.")
    try:
        resolve_view_profile(result.run_directory, checkpoint=result.checkpoint_id)
    except CootViewError as exc:
        print(f"Model/map inspection unavailable for this result: {exc}")
        print("See the full log and report above.")
    else:
        print("Inspect: " + shlex.join([
            "./nasolve", "show", str(result.run_directory),
            "--checkpoint", result.checkpoint_id,
        ]))
    if result.exit_code == 0:
        try:
            run_report = json.loads((run / "report.json").read_text(encoding="utf-8"))
            backbone_policy = requested_backbone_policy(run_report)
            pending_backbone = tuple(backbone_policy["experimental_passthrough_sites"])
        except (OSError, ValueError, BackboneError):
            pending_backbone = ()
        if pending_backbone:
            print()
            print(_color("NON-STANDARD BACKBONE CHEMISTRY — REVIEW REQUIRED", "33"))
            print("Flagged site(s): " + ", ".join(pending_backbone))
            if sys.stdin.isatty():
                # Review is deliberately advisory to refinement status: declining inspection
                # leaves a persistent warning but does not turn a successful refinement into
                # a failed crystallographic run.
                _backbone_review(argparse.Namespace(run=run, coot=None))
            else:
                print("Interactive Coot review was not possible in this session. Run:")
                print("  " + shlex.join(["./nasolve", "backbone-review", str(run)]))
    return result.exit_code


def _refine_doctor(args: argparse.Namespace) -> int:
    try:
        config = load_config()
        run = _workspace_run(args.run, config)
        phenix = discover_phenix(config, explicit=args.phenix_root)
        if args.phenix_root is None:
            remember_phenix(config, phenix)
            save_config(config)
        mtz_dump = phenix.executables.get("phenix.mtz.dump")
        if mtz_dump is None:
            raise RefineDoctorError("Refine Doctor requires phenix.mtz.dump")

        def progress(recipe: str, checkpoint: str, log_path: Path) -> None:
            short_recipe = recipe.rsplit("/", 1)[-1]
            print(
                _color(f"Refine Doctor: {short_recipe} -> {checkpoint}", "36")
                + f"; full output: {log_path}"
            )

        result = execute_refine_doctor(
            run,
            phenix.executables["phenix.refine"],
            mtz_dump,
            phenix_version=phenix.version,
            environment=phenix.environment,
            from_checkpoint=args.from_checkpoint,
            macro_cycles=args.cycles,
            max_trials=args.max_trials,
            progress=progress,
        )
    except (RefineDoctorError, ConfigError, PhenixDiscoveryError) as exc:
        print(f"Refine Doctor error: {exc}", file=sys.stderr)
        return 2
    colors = {
        "REFINE_DOCTOR_GOOD_ENOUGH": "32",
        "REFINE_DOCTOR_RECOMMEND": "36",
        "REFINE_DOCTOR_REVIEW": "33",
        "REFINE_DOCTOR_FLAG_REPAIR_REQUIRED": "31",
    }
    print(f"\n{_color(result.status, colors.get(result.status, '0'))}: {result.message}")
    print(f"Source checkpoint: {result.source_checkpoint}")
    print(
        "Current checkpoint: preserved"
        if result.current_checkpoint_preserved else "Current checkpoint: CHANGED"
    )
    details = result.audit.details
    free_groups = details.get("free_independent_groups")
    fraction = details.get("free_fraction")
    fraction_text = f"{100 * fraction:.1f}%" if isinstance(fraction, (int, float)) else "—"
    print(
        f"Free-R audit: {result.audit.status}; independent free groups "
        f"{free_groups if isinstance(free_groups, int) else '—'} ({fraction_text})"
    )
    for warning in result.audit.warnings:
        print(f"  - {warning}")

    if result.benchmark:
        print("\nStored anomalous benchmark:")
        for item in result.benchmark:
            print(
                f"  {item.get('site')} {item.get('atom_name')}: "
                f"f'' refined {_format_metric(item.get('refined_f_double_prime'))}; "
                f"calculated {_format_metric(item.get('calculated_f_double_prime'))}; "
                f"occupancy {_format_metric(item.get('model_occupancy'))}; "
                f"B {_format_metric(item.get('b_factor'), 2)} A^2 "
                f"[{item.get('source_checkpoint')}]"
            )

    payload = json.loads(result.report_path.read_text(encoding="utf-8"))
    inspection_id = payload.get("inspection_checkpoint")
    triage = payload.get("triage")
    if isinstance(triage, dict):
        print(f"\nTriage: {triage.get('stop_reason')}; {len(result.trials)}/{triage.get('max_trials')} trials")
        for item in triage.get("recipes", []):
            if isinstance(item, dict) and item.get("state") in {"skipped", "failed"}:
                print(f"  {item['state']}: {item.get('recipe')}: {item.get('reason')}")
    candidates = payload.get("candidates")
    if isinstance(candidates, list) and candidates:
        print("\nCheckpoint    Rwork   Rfree     Gap  Recipe")
        for item in candidates:
            if not isinstance(item, dict):
                continue
            marker = ">" if item.get("checkpoint") == result.recommended_checkpoint else " "
            print(
                f"{marker} {str(item.get('checkpoint')):<11} "
                f"{_format_metric(item.get('r_work')):>6} "
                f"{_format_metric(item.get('r_free')):>7} "
                f"{_format_metric(item.get('r_free_minus_r_work')):>7}  "
                f"{item.get('recipe')}"
            )
    print(f"\nRecommendation: {result.recommendation}")
    if not result.recommended_checkpoint and isinstance(inspection_id, str):
        print("Candidate for inspection only (no automatic recommendation):")
        print(f"  nasolve show {shlex.quote(str(result.run_directory))} --checkpoint {shlex.quote(inspection_id)}")
    if isinstance(triage, dict):
        for action in triage.get("next_actions", []):
            print(f"  Next: {action}")
    if result.recommended_checkpoint:
        inspect_command = shlex.join([
            "nasolve", "show", str(result.run_directory),
            "--checkpoint", result.recommended_checkpoint,
        ])
        select_command = shlex.join([
            "nasolve", "checkpoints", "use", str(result.run_directory),
            result.recommended_checkpoint,
        ])
        print(
            "Select after map/model inspection: "
            + select_command
        )
    print(f"Free-R audit log: {result.audit.log_path}")
    print(f"Report: {result.report_path}")
    if result.recommended_checkpoint:
        if sys.stdin.isatty():
            while True:
                try:
                    answer = input(
                        f"\nUse {result.recommended_checkpoint} as the current checkpoint? "
                        "[y/N/i=inspect]: "
                    ).strip().casefold()
                except (EOFError, KeyboardInterrupt):
                    print()
                    answer = "n"
                if answer in {"i", "inspect"}:
                    view_status = _show(argparse.Namespace(
                        target=result.run_directory, dataset=None, stage=None,
                        checkpoint=result.recommended_checkpoint, coot=None,
                    ))
                    if view_status == 0:
                        print("Inspect the model and maps in Coot, then return here to enter y or n.")
                    else:
                        print("Inspection could not open. Retry with i, or inspect later:")
                        print(f"  {inspect_command}")
                    continue
                if answer in {"y", "yes"}:
                    try:
                        selected = select_checkpoint(
                            result.run_directory, result.recommended_checkpoint
                        )
                    except CheckpointError as exc:
                        print(f"Checkpoint error: {exc}", file=sys.stderr)
                        return 2
                    print(_color(f"Current checkpoint: {selected.checkpoint_id}", "36"))
                    break
                if answer in {"", "n", "no"}:
                    print("Current checkpoint unchanged.")
                    print("Inspect without selecting:")
                    print(f"  {inspect_command}")
                    print("Select later:")
                    print(f"  {select_command}")
                    print("Return to the diagnosed source if needed:")
                    print("  " + shlex.join([
                        "nasolve", "checkpoints", "use", str(result.run_directory),
                        result.source_checkpoint,
                    ]))
                    break
                print("Please enter y, n, or i.")
        else:
            print("Inspect without selecting:")
            print(f"  {inspect_command}")
            print("Select later:")
            print(f"  {select_command}")
    return result.exit_code


def _checkpoints(args: argparse.Namespace) -> int:
    try:
        config = load_config()
        if args.checkpoint_action == "use":
            if args.checkpoint is None:
                explicit_run = None
                checkpoint = args.run_or_checkpoint
            else:
                explicit_run = Path(args.run_or_checkpoint)
                checkpoint = args.checkpoint
        else:
            explicit_run = args.run
            checkpoint = None
        run = _workspace_run(explicit_run, config)
        if args.checkpoint_action == "list":
            records, current, bookmarks = list_checkpoints(run)
            names_by_id: dict[str, list[str]] = {}
            for name, target in bookmarks.items():
                names_by_id.setdefault(target, []).append(name)
            print("  Status         Checkpoint   Parent        Rwork/Rfree    Recipe")
            for record in records:
                marker = _status_text(
                    record.status, current=record.checkpoint_id == current
                )
                values = (
                    f"{_format_metric(record.metrics.get('r_work'))}/"
                    f"{_format_metric(record.metrics.get('r_free'))}"
                )
                aliases = names_by_id.get(record.checkpoint_id, [])
                label = f" [{', '.join(aliases)}]" if aliases else ""
                print(
                    f"{marker} {record.checkpoint_id:<12} "
                    f"{(record.parent or '—'):<13} {values:<14} "
                    f"{record.recipe}{label}"
                )
            return 0
        if args.checkpoint_action == "add":
            record = add_checkpoint(
                run,
                name=args.name,
                model=args.model,
                reflections=args.mtz,
                parent=args.from_checkpoint,
            )
            action = "Imported" if args.model else "Bookmarked"
            print(f"{action}: {args.name} -> {record.checkpoint_id}")
            if args.model:
                print("Status: REVIEW (select explicitly after inspection)")
            return 0
        if args.checkpoint_action == "use":
            assert isinstance(checkpoint, str)
            record = select_checkpoint(run, checkpoint)
            print(_color(f"Current checkpoint: {record.checkpoint_id}", "36"))
            print(f"Status: {record.status}")
            print(f"Model: {record.model}")
            return 0
    except (CheckpointError, ConfigError) as exc:
        print(f"Checkpoint error: {exc}", file=sys.stderr)
        return 2
    return 2


def _show(args: argparse.Namespace) -> int:
    try:
        config = load_config()
        if args.target is None:
            run = _workspace_run(None, config)
        else:
            dataset = args.dataset
            if str(args.target).casefold() == "last" and dataset is None:
                dataset = _workspace_dataset(
                    None, config, current_directory_fallback=False
                )
            run = resolve_run(
                args.target,
                dataset,
                active_dataset=(
                    Path(config.workspace.dataset) if config.workspace.dataset else None
                ),
            )
        coot = discover_coot(config, explicit=args.coot)
        remember_coot(config, coot)
        save_config(config)
        result = launch_coot_view(
            run,
            coot.executable,
            stage=args.stage,
            checkpoint=args.checkpoint,
        )
    except (ConfigError, CootDiscoveryError, CootViewError) as exc:
        print(f"Coot view error: {exc}", file=sys.stderr)
        return 2
    print(f"Opened {result.stage} in run: {result.run_directory}")
    print(f"Model source: {result.source}")
    print(f"Model: {result.model_path}")
    print(f"Map source: {result.map_source}")
    print(f"Map: {result.map_path}")
    print("Anomalous map: Coot will also open ANOM/PHANOM at 3 sigma when present")
    for extra in result.extra_model_paths:
        print(f"Additional model: {extra}")
    print(f"Coot working directory: {result.working_directory}")
    print(f"Coot log: {result.log_path}")
    print(f"PID: {result.pid}")
    _print_metal_merge_hint(result.metal_merge)
    return 0


def _relative_to(path: Path | None, root: Path) -> str:
    if path is None:
        return "none"
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _autofdp_files(files: list[Path]) -> tuple[Path, Path]:
    """The model and the MTZ, given in either order."""
    from .autofdp import MODEL_SUFFIXES, AutoFDPError

    models = [path for path in files if path.suffix.lower() in MODEL_SUFFIXES]
    mtzs = [path for path in files if path.suffix.lower() == ".mtz"]
    if len(models) != 1 or len(mtzs) != 1:
        raise AutoFDPError(
            "Give one model (.pdb) and one reflection file (.mtz), not "
            + " and ".join(path.name for path in files)
        )
    return models[0], mtzs[0]


def _print_autofdp_plan(plan) -> None:
    reflections = plan.reflections
    print("AutoFDP inputs")
    print(f"  Model: {plan.model}")
    print(f"  MTZ: {reflections.mtz}")
    print(
        f"  Observations: {','.join(reflections.observation_labels)}; Free-R "
        f"{reflections.free_r_label} (test value {reflections.free_r_test_value})"
    )
    print(
        f"  Wavelength: {reflections.wavelength:g} A "
        f"({reflections.wavelength_source.get('method')})"
    )
    warning = reflections.wavelength_source.get("warning")
    if warning:
        print(f"  Warning: {warning}")
    for path in plan.restraints:
        print(f"  Restraints: {path}")
    theory = "; ".join(
        f"{element} f' {values['f_prime']:.2f}, f'' {values['f_double_prime']:.2f} "
        f"({values['source']})"
        for element, values in plan.theory.items()
    )
    print(f"  Theory: {theory}")
    print("\nAutoFDP step 1: metal B factors from neighbouring nucleotides")
    for site in plan.b_start.sites:
        print(
            f"  {site.site} {site.atom_name} ({site.element}), occupancy "
            f"{site.occupancy:.2f}: B {site.b_before:.2f} -> {site.b_after:.2f}"
        )
        if site.neighbour_b is not None:
            print(f"    {site.neighbour_b.method}:")
            for item in site.neighbour_b.neighbours:
                print(
                    f"      {item.residue} {item.residue_name}  mean B {item.mean_b:.2f}"
                    f"  ({item.atoms} atoms, closest {item.contact})"
                )
        elif site.b_source == "model":
            print("    kept the model's B (--b-start model)")
        else:
            print("    set from --b-start")
    for warning in plan.b_start.warnings:
        print(f"Warning: {warning}")


def _print_fdp_round(result) -> None:
    """Print AutoFDP step 2: refined f'' of every site against theory."""
    round_result = result.round
    status_color = "32" if result.status == "AUTOFDP_FDP_ROUND_DONE" else "31"
    print(
        f"\nAutoFDP step 2: f''-only round ({_relative_to(round_result.directory, result.run_directory)}, "
        f"{round_result.cycles} cycles; only the anomalous scatterers refined, f' fixed at "
        "theory, f'' from 0)"
    )
    print(
        f"R-work {_format_metric(round_result.r_work)}, "
        f"R-free {_format_metric(round_result.r_free)}"
    )
    tolerance = round_result.tolerance
    print(
        f"\n{'Site':<14}{'Occ':>6}{'B':>9}{'f-prime':>9}{'fdp':>9}"
        f"{'theory':>9}{'deviation':>11}"
    )
    for site in round_result.sites:
        refined = "none" if site.f_double_prime_refined is None else f"{site.f_double_prime_refined:.2f}"
        deviation = "n/a" if site.deviation is None else f"{site.deviation:+.2f}"
        occupancy = "-" if site.occupancy is None else f"{site.occupancy:.2f}"
        b_factor = "-" if site.b_factor is None else f"{site.b_factor:.2f}"
        verdict = (
            "" if site.within_tolerance is None
            else f"  within +/-{tolerance:g}" if site.within_tolerance
            else f"  outside +/-{tolerance:g}"
        )
        label = f"{site.site} {site.atom_name}"
        print(
            f"{label:<14}{occupancy:>6}{b_factor:>9}{site.f_prime:>9.2f}{refined:>9}"
            f"{site.f_double_prime_theory:>9.2f}{deviation:>11}{verdict}"
        )
    sources = sorted({f"{site.element} {site.theory_source}" for site in round_result.sites})
    print("Theory: " + "; ".join(sources))
    others = [site for site in round_result.sites if not site.calibrated]
    if others:
        print(
            "Also refined (metal sites left out with --site; B not set): "
            + ", ".join(f"{site.site} {site.atom_name}" for site in others)
        )
    for warning in round_result.warnings:
        print(f"Warning: {warning}")
    print(f"\n{_color(result.status, status_color)}: {result.message}")
    print(f"Run folder: {result.run_directory}")


def _open_fdp_view(args: argparse.Namespace, result) -> None:
    """Open the f'' round's model and maps; a Coot problem never fails the run."""
    round_result = result.round
    if round_result.model_path is None or round_result.map_coefficients is None:
        return
    if args.no_coot:
        print(f"Model: {round_result.model_path}")
        print(f"Maps: {round_result.map_coefficients}")
        return
    try:
        config = load_config()
        coot = discover_coot(config, explicit=args.coot)
        remember_coot(config, coot)
        save_config(config)
        view = launch_coot_files(
            round_result.model_path,
            round_result.map_coefficients,
            coot.executable,
            working_directory=result.run_directory / "CootGUI",
            source=f"AutoFDP f''-only round ({result.run_directory.name})",
            dictionaries=result.plan.restraints,
        )
    except (ConfigError, CootDiscoveryError, CootViewError) as exc:
        print(f"Coot view error: {exc}", file=sys.stderr)
        print(f"Model: {round_result.model_path}")
        print(f"Maps: {round_result.map_coefficients}")
        return
    print(f"Opened Coot (PID {view.pid}) with {view.model_path.name} and {view.map_path.name}")
    print("Anomalous map: Coot also opens ANOM/PHANOM at 3 sigma when present")
    print(f"Coot log: {view.log_path}")


def _autofdp(args: argparse.Namespace) -> int:
    """AutoFDP steps 1 and 2 from a model and an MTZ."""
    from .autofdp import (
        STATUS_ROUND_DONE,
        AutoFDPError,
        execute_autofdp,
        parse_b_start,
        parse_theory_overrides,
        plan_autofdp,
    )

    if args.element:
        print("AutoFDP error: --element is not implemented yet", file=sys.stderr)
        return 2
    try:
        model, mtz = _autofdp_files([args.model, args.mtz])
        b_start = parse_b_start(args.b_start)
        theory = parse_theory_overrides(args.theory)
        if args.cycles < 1:
            raise AutoFDPError("--cycles must be at least 1")
        if not math.isfinite(args.tolerance) or args.tolerance <= 0:
            raise AutoFDPError("--tolerance must be a positive number of electrons")
        if args.coot is not None and args.no_coot:
            raise AutoFDPError("--coot and --no-coot cannot be combined")
        config = load_config()
        phenix = discover_phenix(config, explicit=args.phenix_root)
        if args.phenix_root is None:
            remember_phenix(config, phenix)
            save_config(config)
        refine = phenix.executables.get("phenix.refine")
        mtz_dump = phenix.executables.get("phenix.mtz.dump")
        if refine is None or mtz_dump is None:
            raise AutoFDPError("AutoFDP requires phenix.refine and phenix.mtz.dump")
        plan = plan_autofdp(
            model, mtz, refine, mtz_dump,
            sites=tuple(args.site), b_start=b_start, theory_overrides=theory,
            wavelength=args.wavelength, restraints=tuple(args.cif),
            output_root=args.output, environment=phenix.environment,
        )
        _print_autofdp_plan(plan)
        if args.plan:
            print(
                f"\nThen: one f''-only phenix.refine round of {args.cycles} cycles (every "
                "strategy off except the anomalous scatterers; f' fixed at theory, "
                f"f'' from 0), in a new run folder of {plan.output_root}"
            )
            print("Plan only: no refinement was run and nothing was written.")
            return 0
        print(f"\nAutoFDP f'' round started ({args.cycles} cycles)")
        result = execute_autofdp(
            plan, refine,
            phenix_version=phenix.version, b_start=b_start, sites=tuple(args.site),
            cycles=args.cycles, tolerance=args.tolerance, environment=phenix.environment,
        )
    except (AutoFDPError, ConfigError, PhenixDiscoveryError) as exc:
        print(f"AutoFDP error: {exc}", file=sys.stderr)
        return 1
    _print_fdp_round(result)
    _open_fdp_view(args, result)
    print("Occupancy adjustment and the later AutoFDP steps are not implemented yet.")
    return 0 if result.status == STATUS_ROUND_DONE else 1


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "configure" and args.configure_target == "phenix":
        return _configure_phenix(args.path)
    if args.command == "configure" and args.configure_target == "coot":
        return _configure_coot(args.path)
    if args.command == "workspace":
        return _workspace(args)
    if args.command == "preset":
        return _preset(args)
    if args.command == "campaign":
        return _campaign(args)
    if args.command == "check":
        return _check(args.phenix_root)
    if args.command == "automr":
        return _automr(args)
    if args.command == "postmr":
        return _postmr(args)
    if args.command == "backbone-review":
        return _backbone_review(args)
    if args.command == "autosol":
        return _autosol(args)
    if args.command == "autorefine":
        return _autorefine(args)
    if args.command in {"refine-doctor", "refine_doctor"}:
        return _refine_doctor(args)
    if args.command == "autofdp":
        return _autofdp(args)
    if args.command == "checkpoints":
        return _checkpoints(args)
    if args.command == "show":
        return _show(args)
    return 2
