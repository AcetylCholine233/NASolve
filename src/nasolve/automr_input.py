"""Dataset discovery and shared standard/nonstandard AutoMR input schema."""

from __future__ import annotations

import configparser
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Collection, Mapping

from .phosphate import PhosphateError, validate_op3_sites, validate_phosphate_intent
from .backbone import BackboneError, validate_backbone_sites
from .presets import PresetError, ProjectPreset, load_preset
from .reflection_data import DEFAULT_PHASER_DATA, ReflectionDataError, normalize_phaser_data
from .residue_aliases import LigandCodeError, ResolvedLigand, resolve_ligand, resolve_pair


class AutoMRInputError(RuntimeError):
    """Raised when an AutoMR request is incomplete, unsafe, or ambiguous."""


@dataclass(frozen=True)
class FrameSpec:
    name: str
    accession: str
    aliases: tuple[str, ...]
    catalogue_directory: str
    fallback_model: str
    recipe_preset: str | None = None


FRAME_SPECS: tuple[FrameSpec, ...] = (
    FrameSpec("W", "5W6W", ("W", "5W6W"), "5W6W", "C_G.pdb", "5w6w"),
    FrameSpec("3GBI", "3GBI", ("3GBI",), "3GBI", "C_C.pdb"),
)


@dataclass
class AutoMRIntent:
    mode: str | None = None
    frame: str | None = None
    pair: str | None = None
    model: str | None = None
    model_family: str | None = None
    sequence_file: str | None = None
    mirror: bool = False
    allow_p1_standard: bool = False
    sequences: dict[str, str] = field(default_factory=dict)
    mutations: dict[str, str] = field(default_factory=dict)
    source: Path | None = None
    allow_op3_sites: tuple[str, ...] = ()
    op3_sites_explicit: bool = False
    backbone_sites: dict[str, str] = field(default_factory=dict)
    allow_unreviewed_backbone: bool = False
    sequence_reference: str | None = None
    phaser_data: str | None = None


@dataclass(frozen=True)
class DatasetFiles:
    root: Path
    reflections: Path
    metadata: Path
    summary: Path


@dataclass(frozen=True)
class ResolvedAutoMRInput:
    dataset: DatasetFiles
    mode: str
    frame: FrameSpec | None
    pair_text: str | None
    pair: tuple[ResolvedLigand, ResolvedLigand] | None
    model: Path
    model_source: str
    model_pair: tuple[ResolvedLigand, ResolvedLigand] | None
    exact_pair_model: bool | None
    catalogue_warnings: tuple[str, ...]
    allow_p1_standard: bool
    mirror: bool
    sequences: dict[str, str]
    sequence_file: Path | None
    mutations: dict[str, ResolvedLigand]
    config_source: Path | None
    allow_op3_sites: tuple[str, ...] = ()
    phosphate_intent: dict[str, object] | None = None
    backbone_sites: dict[str, str] = field(default_factory=dict)
    allow_unreviewed_backbone: bool = False
    sequence_reference: Path | None = None
    sequence_reference_label: str | None = None
    sequence_thread: dict[str, object] | None = None
    model_selector: str | None = None
    model_family: str | None = None
    model_provider: dict[str, object] | None = None
    frame_sequence_source: Path | None = None
    phaser_data: str = DEFAULT_PHASER_DATA


_ALLOWED_SECTIONS = {"automr", "sequences", "mutations", "backbones"}
_ALLOWED_AUTOMR_KEYS = {
    "mode",
    "frame",
    "pair",
    "model",
    "model_family",
    "sequence_file",
    "sequence_reference",
    "mirror",
    "allow_p1_standard",
    "allow_op3_sites",
    "five_prime_phosphate_sites",
    "allow_unreviewed_backbone",
    "phaser_data",
}


def _validated_model_family(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise AutoMRInputError("model_family must be a stable text identifier")
    family = value.strip()
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", family) is None:
        raise AutoMRInputError(
            "model_family must use 1-64 letters, numbers, dot, underscore or hyphen"
        )
    return family


def _validated_sequences(
    sequences: Mapping[str, str],
    context: str,
) -> dict[str, str]:
    validated: dict[str, str] = {}
    for raw_chain, raw_sequence in sequences.items():
        chain = raw_chain.strip()
        sequence = "".join(raw_sequence.split()).upper()
        if not chain or not sequence:
            raise AutoMRInputError(f"{context} chain names and sequences cannot be empty")
        if len(chain) != 1:
            raise AutoMRInputError(
                f"{context} chain {chain!r} must be one PDB chain identifier"
            )
        invalid = sorted(set(sequence) - set("ACGTU"))
        if invalid:
            raise AutoMRInputError(
                f"{context} chain {chain} contains unsupported sequence symbol(s): "
                + ", ".join(invalid)
            )
        if chain in validated:
            raise AutoMRInputError(f"{context} contains duplicate chain {chain!r}")
        validated[chain] = sequence
    return validated


def read_sequence_file(path: Path) -> dict[str, str]:
    """Read chain-labelled FASTA or ``CHAIN = SEQUENCE`` text."""
    source = path.expanduser().resolve()
    try:
        lines = source.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise AutoMRInputError(f"Could not read sequence file {source}: {exc}") from exc
    meaningful = [line.strip() for line in lines if line.strip() and not line.lstrip().startswith("#")]
    if not meaningful:
        raise AutoMRInputError(f"Sequence file is empty: {source}")
    parsed: dict[str, str] = {}
    if any(line.startswith(">") for line in meaningful):
        current: str | None = None
        chunks: list[str] = []

        def finish() -> None:
            if current is not None:
                if current in parsed:
                    raise AutoMRInputError(
                        f"Sequence file contains duplicate FASTA chain {current!r}"
                    )
                parsed[current] = "".join(chunks)

        for line in meaningful:
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                current = header.split()[0] if header else None
                chunks = []
                if current is None:
                    raise AutoMRInputError("FASTA sequence header must name a PDB chain")
            else:
                if current is None:
                    raise AutoMRInputError(
                        "FASTA sequence data appeared before the first chain header"
                    )
                chunks.append(line)
        finish()
    else:
        for line in meaningful:
            if line.count("=") != 1:
                raise AutoMRInputError(
                    "Non-FASTA sequence files must use CHAIN = SEQUENCE on each line"
                )
            chain, sequence = line.split("=", 1)
            chain = chain.strip()
            if chain in parsed:
                raise AutoMRInputError(
                    f"Sequence file contains duplicate chain {chain!r}"
                )
            parsed[chain] = sequence.strip()
    return _validated_sequences(parsed, f"Sequence file {source.name}")


def _parser() -> configparser.ConfigParser:
    parser = configparser.ConfigParser(
        interpolation=None,
        delimiters=("=",),
        comment_prefixes=("#", ";"),
        inline_comment_prefixes=("#", ";"),
        strict=True,
        empty_lines_in_values=False,
    )
    parser.optionxform = str
    return parser


def read_intent(path: Path | None) -> AutoMRIntent:
    """Read the canonical INI-style ``nasolve.txt`` input file."""
    if path is None:
        return AutoMRIntent()
    source = path.expanduser().resolve()
    if not source.is_file():
        raise AutoMRInputError(f"NASolve input file does not exist: {source}")
    parser = _parser()
    try:
        with source.open(encoding="utf-8") as handle:
            parser.read_file(handle)
    except (OSError, configparser.Error) as exc:
        raise AutoMRInputError(f"Could not read {source}: {exc}") from exc
    unknown_sections = set(parser.sections()) - _ALLOWED_SECTIONS
    if unknown_sections:
        raise AutoMRInputError(
            "Unknown nasolve.txt section(s): " + ", ".join(sorted(unknown_sections))
        )
    if "automr" not in parser:
        raise AutoMRInputError("nasolve.txt must contain an [automr] section")
    automr = dict(parser["automr"])
    unknown_keys = set(automr) - _ALLOWED_AUTOMR_KEYS
    if unknown_keys:
        raise AutoMRInputError(
            "Unknown [automr] setting(s): " + ", ".join(sorted(unknown_keys))
        )
    try:
        allow_p1_standard = parser["automr"].getboolean(
            "allow_p1_standard", fallback=False
        )
    except ValueError as exc:
        raise AutoMRInputError(
            "[automr] allow_p1_standard must be true or false"
        ) from exc
    try:
        mirror = parser["automr"].getboolean("mirror", fallback=False)
    except ValueError as exc:
        raise AutoMRInputError("[automr] mirror must be true or false") from exc
    try:
        allow_unreviewed_backbone = parser["automr"].getboolean(
            "allow_unreviewed_backbone", fallback=False
        )
    except ValueError as exc:
        raise AutoMRInputError(
            "[automr] allow_unreviewed_backbone must be true or false"
        ) from exc
    sequences = _validated_sequences({
        chain.strip(): "".join(sequence.split())
        for chain, sequence in parser["sequences"].items()
    } if "sequences" in parser else {}, "[sequences]")
    mutations = {
        site.strip(): residue.strip()
        for site, residue in parser["mutations"].items()
    } if "mutations" in parser else {}
    if any(not site or not residue for site, residue in mutations.items()):
        raise AutoMRInputError("Mutation sites and residue values cannot be empty")
    if "allow_op3_sites" in automr and "five_prime_phosphate_sites" in automr:
        raise AutoMRInputError(
            "Use five_prime_phosphate_sites (preferred) or legacy allow_op3_sites, not both"
        )
    phosphate_key = (
        "five_prime_phosphate_sites"
        if "five_prime_phosphate_sites" in automr
        else "allow_op3_sites"
    )
    raw_op3 = automr.get(phosphate_key, "").strip()
    try:
        allow_op3 = validate_op3_sites(
            [v.strip() for v in raw_op3.split(",")] if raw_op3 else []
        )
        backbone_sites = validate_backbone_sites(
            {site.strip(): mode.strip() for site, mode in parser["backbones"].items()}
            if "backbones" in parser else {}
        )
    except (PhosphateError, BackboneError) as exc:
        raise AutoMRInputError(str(exc)) from exc
    return AutoMRIntent(
        mode=automr.get("mode") or None,
        frame=automr.get("frame") or None,
        pair=automr.get("pair") or None,
        model=automr.get("model") or None,
        model_family=_validated_model_family(automr.get("model_family")),
        sequence_file=automr.get("sequence_file") or None,
        sequence_reference=automr.get("sequence_reference"),
        mirror=mirror,
        allow_p1_standard=allow_p1_standard,
        sequences=sequences,
        mutations=mutations,
        source=source,
        allow_op3_sites=allow_op3,
        op3_sites_explicit=(
            "allow_op3_sites" in automr or "five_prime_phosphate_sites" in automr
        ),
        backbone_sites=backbone_sites,
        allow_unreviewed_backbone=allow_unreviewed_backbone,
        phaser_data=_validated_phaser_data(automr.get("phaser_data"), "[automr] phaser_data"),
    )


def _validated_phaser_data(value: object, label: str) -> str | None:
    if value is None:
        return None
    try:
        return normalize_phaser_data(value)
    except ReflectionDataError as exc:
        raise AutoMRInputError(f"{label}: {exc}") from exc


def discover_dataset(root: Path) -> DatasetFiles:
    dataset = root.expanduser().resolve()
    if not dataset.is_dir():
        raise AutoMRInputError(f"Dataset directory does not exist: {dataset}")

    mtz_files = sorted(
        path for path in dataset.iterdir()
        if path.is_file() and path.suffix.casefold() == ".mtz"
    )
    staraniso_candidates = [
        path for path in mtz_files
        if "staraniso" in _normalized_name(path) and "alldata" in _normalized_name(path)
    ]
    if len(staraniso_candidates) == 1:
        reflections = staraniso_candidates[0]
    elif len(staraniso_candidates) > 1:
        raise AutoMRInputError(
            "Ambiguous STARANISO all-data MTZ files: "
            + ", ".join(path.name for path in staraniso_candidates)
        )
    elif len(mtz_files) == 1:
        reflections = mtz_files[0]
    elif not mtz_files:
        raise AutoMRInputError("Missing required reflection MTZ file")
    else:
        raise AutoMRInputError(
            "Ambiguous reflection MTZ files: " + ", ".join(path.name for path in mtz_files)
        )

    metadata_candidates = sorted(
        path for path in dataset.iterdir()
        if path.is_file()
        and path.suffix.casefold() == ".cif"
        and _normalized_name(path).startswith("data1")
    )
    if not metadata_candidates:
        raise AutoMRInputError("Missing required Data_1*.cif metadata file")
    if len(metadata_candidates) > 1:
        raise AutoMRInputError(
            "Ambiguous Data_1*.cif metadata files: "
            + ", ".join(path.name for path in metadata_candidates)
        )
    summaries = sorted(
        path for path in dataset.iterdir()
        if path.is_file() and path.name.casefold() == "summary.html"
    )
    if not summaries:
        raise AutoMRInputError("Missing required file: summary.html")
    if len(summaries) > 1:
        raise AutoMRInputError(
            "Ambiguous summary files: " + ", ".join(path.name for path in summaries)
        )
    summary = summaries[0]
    return DatasetFiles(dataset, reflections, metadata_candidates[0], summary)


def _normalized_name(path: Path) -> str:
    """Return a case/punctuation-insensitive filename stem for discovery."""
    return re.sub(r"[^a-z0-9]+", "", path.stem.casefold())


def normalize_frame(value: str) -> FrameSpec:
    cleaned = value.strip().casefold()
    for spec in FRAME_SPECS:
        if cleaned in {alias.casefold() for alias in spec.aliases}:
            return spec
    supported = ", ".join(spec.name for spec in FRAME_SPECS)
    raise AutoMRInputError(f"Unknown standard frame {value!r}; choose {supported}")


def locate_frames_directory(
    explicit: Path | None = None,
    environ: Mapping[str, str] | None = None,
) -> Path:
    """Find the development or future package location of approved MR frames."""
    env = os.environ if environ is None else environ
    candidates: list[tuple[str, Path]] = []
    if explicit is not None:
        candidates.append(("--frames-dir", explicit.expanduser()))
    if env.get("NASOLVE_MR_FRAMES"):
        candidates.append(("NASOLVE_MR_FRAMES", Path(env["NASOLVE_MR_FRAMES"]).expanduser()))
    package_dir = Path(__file__).resolve().parent
    candidates.extend([
        ("package data", package_dir / "data" / "MR_frames"),
        ("repository MR_frames", package_dir.parents[1] / "MR_frames"),
    ])
    for _, candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    attempted = ", ".join(f"{label}={path}" for label, path in candidates)
    raise AutoMRInputError(
        "Standard mode needs the approved MR_frames directory. "
        "Use --frames-dir PATH or set NASOLVE_MR_FRAMES. Tried: " + attempted
    )


def _pair_from_model_name(
    path: Path,
    valid_ligand_codes: Collection[str] | None,
) -> tuple[ResolvedLigand, ResolvedLigand]:
    """Resolve the ordered TOKEN_TOKEN pair encoded by a catalogue filename."""
    if path.stem.count("_") != 1:
        raise AutoMRInputError(
            f"Invalid standard model filename {path.name!r}; use FIRST_SECOND.pdb"
        )
    first, second = path.stem.split("_", 1)
    try:
        return (
            resolve_ligand(first, valid_ligand_codes),
            resolve_ligand(second, valid_ligand_codes),
        )
    except LigandCodeError as exc:
        raise AutoMRInputError(f"Invalid standard model filename {path.name!r}: {exc}") from exc


def _resolve_frame_model(
    spec: FrameSpec,
    frames_dir: Path,
    requested_pair: tuple[ResolvedLigand, ResolvedLigand],
    valid_ligand_codes: Collection[str] | None,
) -> tuple[Path, tuple[ResolvedLigand, ResolvedLigand], bool, tuple[str, ...]]:
    catalogue = frames_dir / spec.catalogue_directory
    if not catalogue.is_dir():
        raise AutoMRInputError(
            f"Missing standard frame catalogue: {catalogue}"
        )
    files = sorted(
        path for path in catalogue.iterdir()
        if path.is_file() and path.suffix.casefold() == ".pdb"
    )
    if not files:
        raise AutoMRInputError(f"No PDB models found in standard frame catalogue {catalogue}")

    requested_codes = tuple(ligand.ligand_code for ligand in requested_pair)
    indexed: list[tuple[Path, tuple[ResolvedLigand, ResolvedLigand]]] = []
    warnings: list[str] = []
    for path in files:
        try:
            indexed.append((path, _pair_from_model_name(path, valid_ligand_codes)))
        except AutoMRInputError as exc:
            warnings.append(str(exc))
    matches = [
        (path, pair)
        for path, pair in indexed
        if tuple(ligand.ligand_code for ligand in pair) == requested_codes
    ]
    if len(matches) > 1:
        raise AutoMRInputError(
            "Ambiguous exact standard models for "
            f"{requested_codes[0]}:{requested_codes[1]}: "
            + ", ".join(path.name for path, _ in matches)
        )
    if matches:
        path, pair = matches[0]
        return path.resolve(), pair, True, tuple(warnings)

    fallback = catalogue / spec.fallback_model
    if not fallback.is_file():
        raise AutoMRInputError(
            f"Required fallback model for frame {spec.name} is missing: {fallback}"
        )
    fallback_pair = _pair_from_model_name(fallback, valid_ligand_codes)
    return fallback.resolve(), fallback_pair, False, tuple(warnings)


def _resolve_standard_model_override(
    dataset: Path,
    spec: FrameSpec,
    frames_dir: Path,
    configured: str,
) -> tuple[Path, dict[str, object]]:
    """Resolve one explicit standard-frame model without inferring compatibility."""
    selector = configured.strip()
    if not selector or any(char in selector for char in "\r\n\x00"):
        raise AutoMRInputError("Explicit standard model selector is malformed")
    relative = Path(selector)
    if (
        relative.is_absolute()
        or "\\" in selector
        or any(part in {"", ".", ".."} for part in relative.parts)
        or relative.suffix.casefold() != ".pdb"
    ):
        raise AutoMRInputError(
            "Standard model override must be a safe relative PDB path"
        )

    catalogue = (frames_dir / spec.catalogue_directory).resolve()
    if not catalogue.is_dir():
        raise AutoMRInputError(f"Missing standard frame catalogue: {catalogue}")

    candidates: list[tuple[str, Path]] = []
    dataset_model = (dataset / relative).resolve()
    try:
        dataset_model.relative_to(dataset)
    except ValueError as exc:
        raise AutoMRInputError("Standard model override must remain inside the dataset or frame catalogue") from exc
    if dataset_model.is_file():
        candidates.append(("dataset", dataset_model))

    frame_model = (catalogue / relative).resolve()
    try:
        frame_model.relative_to(catalogue)
    except ValueError as exc:
        raise AutoMRInputError("Standard model override must remain inside the dataset or frame catalogue") from exc
    if frame_model.is_file():
        candidates.append(("frame-catalogue", frame_model))

    if not candidates:
        raise AutoMRInputError(
            f"Explicit standard model {selector!r} was not found in the dataset "
            f"or {catalogue}"
        )
    if len(candidates) > 1:
        raise AutoMRInputError(
            f"Ambiguous explicit standard model {selector!r}: present in both "
            "the dataset and frame catalogue"
        )
    location, model = candidates[0]
    return model, {
        "kind": "explicit-standard-model",
        "selection": "user-forced",
        "frame": spec.name,
        "location": location,
        "selector": selector,
    }


def _resolve_dataset_model(dataset: Path, configured: str | None) -> tuple[Path, str]:
    if configured:
        relative = Path(configured).expanduser()
        if relative.is_absolute():
            raise AutoMRInputError("Nonstandard model paths in nasolve.txt must be dataset-relative")
        model = (dataset / relative).resolve()
        try:
            model.relative_to(dataset)
        except ValueError as exc:
            raise AutoMRInputError("Nonstandard model path must remain inside the dataset") from exc
        if not model.is_file():
            raise AutoMRInputError(f"Configured MR model does not exist: {model}")
        if model.suffix.casefold() != ".pdb":
            raise AutoMRInputError("The current nonstandard model layer accepts PDB files only")
        return model, "nasolve.txt"
    candidates = sorted(path.resolve() for path in dataset.glob("*.pdb") if path.is_file())
    if not candidates:
        raise AutoMRInputError(
            "MODEL_REQUIRED: no top-level PDB model was found and standard mode was not selected"
        )
    if len(candidates) > 1:
        raise AutoMRInputError(
            "Ambiguous nonstandard MR models: " + ", ".join(path.name for path in candidates)
        )
    return candidates[0], "dataset discovery"


def _resolve_relative_input(dataset: Path, configured: str, label: str) -> Path:
    relative = Path(configured).expanduser()
    if relative.is_absolute():
        raise AutoMRInputError(f"{label} paths in nasolve.txt must be dataset-relative")
    selected = (dataset / relative).resolve()
    try:
        selected.relative_to(dataset)
    except ValueError as exc:
        raise AutoMRInputError(f"{label} path must remain inside the dataset") from exc
    if not selected.is_file():
        raise AutoMRInputError(f"Configured {label.lower()} does not exist: {selected}")
    return selected


def _validate_mutation_site(site: str) -> None:
    if site.count(":") != 1:
        raise AutoMRInputError(
            f"Mutation site {site!r} must use CHAIN:RESID syntax, for example A:8"
        )
    chain, residue = (part.strip() for part in site.split(":", 1))
    if not chain or not residue:
        raise AutoMRInputError(
            f"Mutation site {site!r} must use CHAIN:RESID syntax, for example A:8"
        )


def resolve_automr_input(
    root: Path,
    intent: AutoMRIntent,
    frame_override: str | None = None,
    pair_override: str | None = None,
    frames_dir: Path | None = None,
    allow_p1_standard: bool = False,
    mirror_override: bool = False,
    model_override: str | None = None,
    environ: Mapping[str, str] | None = None,
    valid_ligand_codes: Collection[str] | None = None,
    *,
    recipe: ProjectPreset | None = None,
    phaser_data_override: str | None = None,
) -> ResolvedAutoMRInput:
    """Apply CLI precedence, validate intent, and select exactly one model."""
    try:
        allow_op3 = validate_op3_sites(intent.allow_op3_sites)
        backbone_sites = validate_backbone_sites(intent.backbone_sites)
    except (PhosphateError, BackboneError) as exc:
        raise AutoMRInputError(str(exc)) from exc
    if type(intent.allow_unreviewed_backbone) is not bool:
        raise AutoMRInputError("allow_unreviewed_backbone must be true or false")
    dataset = discover_dataset(root)
    configured_mode = intent.mode.strip().casefold() if intent.mode else None
    if configured_mode not in {None, "standard", "nonstandard"}:
        raise AutoMRInputError("[automr] mode must be standard or nonstandard")

    selected_frame = frame_override or intent.frame
    # A command-line frame explicitly changes the effective mode.  Otherwise an
    # explicit nonstandard mode is not silently reinterpreted.
    if frame_override:
        mode = "standard"
    elif configured_mode:
        mode = configured_mode
    elif selected_frame:
        mode = "standard"
    else:
        mode = "nonstandard"

    pair_text = pair_override or intent.pair
    frame: FrameSpec | None = None
    pair: tuple[ResolvedLigand, ResolvedLigand] | None = None
    model_pair: tuple[ResolvedLigand, ResolvedLigand] | None = None
    exact_pair_model: bool | None = None
    catalogue_warnings: tuple[str, ...] = ()
    model_selector = model_override if model_override is not None else intent.model
    model_family = _validated_model_family(intent.model_family)
    if model_family is not None:
        if not isinstance(intent.model, str) or not intent.model.strip():
            raise AutoMRInputError(
                "[automr] model_family requires an explicit model = selector"
            )
        if model_override is not None and model_override != intent.model:
            raise AutoMRInputError(
                "Configured model_family is bound to model = "
                f"{intent.model!r} and cannot follow a different --model override"
            )
    model_provider: dict[str, object] | None = None
    frame_sequence_source: Path | None = None
    effective_allow_p1 = bool(allow_p1_standard or intent.allow_p1_standard)
    # The command line overrides nasolve.txt; neither means F (amplitudes).
    phaser_data = (
        _validated_phaser_data(phaser_data_override, "--phaser-data")
        if phaser_data_override is not None
        else _validated_phaser_data(intent.phaser_data, "[automr] phaser_data")
    ) or DEFAULT_PHASER_DATA
    effective_mirror = bool(mirror_override or intent.mirror)
    if mode == "standard":
        if not selected_frame:
            raise AutoMRInputError("Standard mode requires frame = W or frame = 3GBI")
        if not pair_text:
            raise AutoMRInputError("Standard mode requires an ordered pair, for example pair = D:T")
        frame = normalize_frame(selected_frame)
        try:
            pair = resolve_pair(pair_text, valid_ligand_codes)
        except LigandCodeError as exc:
            raise AutoMRInputError(str(exc)) from exc
        located_frames = locate_frames_directory(frames_dir, environ)
        catalogue = located_frames / frame.catalogue_directory
        sequence = catalogue / "seq_base.txt"
        frame_sequence_source = sequence.resolve() if sequence.is_file() else None
        if model_selector is not None:
            model, model_provider = _resolve_standard_model_override(
                dataset.root, frame, located_frames, model_selector
            )
            model_pair = None
            exact_pair_model = False
            model_source = (
                f"explicit standard model ({frame.name}; "
                f"{model_provider['location']}:{model_provider['selector']})"
            )
        else:
            model, model_pair, exact_pair_model, catalogue_warnings = _resolve_frame_model(
                frame, located_frames, pair, valid_ligand_codes
            )
            model_provider = {
                "kind": "standard-frame-catalogue",
                "selection": "exact-pair" if exact_pair_model else "fallback",
                "frame": frame.name,
                "location": "frame-catalogue",
                "selector": model.name,
            }
            if exact_pair_model:
                model_source = f"standard frame catalogue ({frame.name}; exact pair)"
            else:
                model_source = (
                    f"standard frame catalogue ({frame.name}; fallback {frame.fallback_model})"
                )
    else:
        if selected_frame and not frame_override:
            raise AutoMRInputError("Nonstandard mode cannot also specify a standard frame")
        if pair_text:
            raise AutoMRInputError(
                "pair = is only valid in standard mode; use exact [mutations] sites for nonstandard models"
            )
        model, model_source = _resolve_dataset_model(dataset.root, model_selector)
        model_provider = {
            "kind": (
                "explicit-nonstandard-model" if model_selector is not None
                else "discovered-nonstandard-model"
            ),
            "selection": "user-forced" if model_selector is not None else "single-pdb-discovery",
            "location": "dataset",
            "selector": model_selector or model.relative_to(dataset.root).as_posix(),
        }

    if model_family is not None:
        assert model_provider is not None
        model_provider["construct_family"] = model_family

    sequence_file: Path | None = None
    sequences = dict(intent.sequences)
    if intent.sequence_file:
        if mode != "nonstandard":
            raise AutoMRInputError("sequence_file is currently supported only in nonstandard mode")
        if sequences:
            raise AutoMRInputError(
                "Use either [automr] sequence_file or inline [sequences], not both"
            )
        sequence_file = _resolve_relative_input(
            dataset.root, intent.sequence_file, "Sequence file"
        )
        sequences = read_sequence_file(sequence_file)

    sequence_reference: Path | None = None
    if intent.sequence_reference is not None:
        requested_reference = intent.sequence_reference
        if not isinstance(requested_reference, str) or not requested_reference.strip():
            raise AutoMRInputError("sequence_reference must be an explicit reference ID or dataset-relative file")
        requested_reference = requested_reference.strip()
        if requested_reference == "w-metal-scaffold":
            sequence_reference = (
                Path(__file__).resolve().parent / "data" / "sequence_references"
                / "w-metal-scaffold.json"
            )
            if not sequence_reference.is_file():
                raise AutoMRInputError("The selected w-metal-scaffold sequence reference is missing")
        else:
            sequence_reference = _resolve_relative_input(
                dataset.root, requested_reference, "Sequence reference"
            )

    resolved_mutations: dict[str, ResolvedLigand] = {}
    for site, residue in intent.mutations.items():
        _validate_mutation_site(site)
        try:
            resolved_mutations[site] = resolve_ligand(residue, valid_ligand_codes)
        except LigandCodeError as exc:
            raise AutoMRInputError(f"Mutation {site}: {exc}") from exc

    # A chosen recipe is itself an explicit chemistry request, not model inference.
    # Campaign planning passes its selected card. Execution uses frozen values and
    # deliberately never calls this resolver again.
    declaration = None
    if frame is not None:
        try:
            selected_recipe = recipe or (load_preset(frame.recipe_preset) if frame.recipe_preset else None)
        except PresetError as exc:
            raise AutoMRInputError(f"Cannot load selected frame recipe: {exc}") from exc
        if selected_recipe is not None:
            if normalize_frame(selected_recipe.automr_defaults["frame"]).name != frame.name:
                raise AutoMRInputError("Selected recipe and effective frame disagree")
            declaration = selected_recipe.phosphate_declaration()
    elif recipe is not None:
        raise AutoMRInputError("A standard recipe cannot authorize a nonstandard model")
    explicit = bool(intent.op3_sites_explicit or intent.allow_op3_sites)
    if not explicit and declaration is not None:
        allow_op3 = tuple(declaration["terminal_phosphate_sites"])
    phosphate_intent = {
        "schema_version": 1,
        "source": "dataset" if explicit else ("recipe" if declaration is not None else "none"),
        "recipe": declaration,
        "allow_op3_sites": list(allow_op3),
    }
    try:
        validate_phosphate_intent(phosphate_intent, allow_op3)
    except PhosphateError as exc:
        raise AutoMRInputError(str(exc)) from exc

    return ResolvedAutoMRInput(
        dataset=dataset,
        mode=mode,
        frame=frame,
        pair_text=pair_text,
        pair=pair,
        model=model,
        model_source=model_source,
        model_pair=model_pair,
        exact_pair_model=exact_pair_model,
        catalogue_warnings=catalogue_warnings,
        allow_p1_standard=effective_allow_p1,
        mirror=effective_mirror,
        sequences=sequences,
        sequence_file=sequence_file,
        sequence_reference=sequence_reference,
        sequence_reference_label=requested_reference if sequence_reference is not None else None,
        model_selector=model_selector,
        model_family=model_family,
        model_provider=model_provider,
        frame_sequence_source=frame_sequence_source,
        mutations=resolved_mutations,
        config_source=intent.source,
        allow_op3_sites=allow_op3,
        phosphate_intent=phosphate_intent,
        backbone_sites=backbone_sites,
        allow_unreviewed_backbone=intent.allow_unreviewed_backbone,
        phaser_data=phaser_data,
    )


def format_intent(resolved: ResolvedAutoMRInput) -> str:
    """Return a canonical, re-readable snapshot of effective user intent."""
    lines = ["[automr]", f"mode = {resolved.mode}"]
    if resolved.mode == "standard":
        assert resolved.frame is not None and resolved.pair_text is not None
        lines.extend([f"frame = {resolved.frame.name}", f"pair = {resolved.pair_text}"])
        if resolved.model_selector is not None:
            lines.append(f"model = {resolved.model_selector}")
        if resolved.model_family is not None:
            lines.append(f"model_family = {resolved.model_family}")
        if resolved.allow_p1_standard:
            lines.append("allow_p1_standard = true")
    else:
        relative_model = (
            resolved.model_selector
            if resolved.model_selector is not None
            else resolved.model.relative_to(resolved.dataset.root).as_posix()
        )
        lines.append(f"model = {relative_model}")
        if resolved.model_family is not None:
            lines.append(f"model_family = {resolved.model_family}")
    if resolved.sequence_reference is not None:
        if resolved.sequence_reference_label is not None:
            reference_text = resolved.sequence_reference_label.strip()
            if not reference_text or any(char in reference_text for char in "\r\n\x00"):
                raise AutoMRInputError("Sequence reference label is malformed")
        else:
            builtin = Path(__file__).resolve().parent / "data" / "sequence_references" / "w-metal-scaffold.json"
            selected = resolved.sequence_reference.resolve()
            if selected == builtin.resolve():
                reference_text = "w-metal-scaffold"
            else:
                try:
                    reference_text = selected.relative_to(resolved.dataset.root.resolve()).as_posix()
                except ValueError as exc:
                    raise AutoMRInputError("Sequence reference must be built-in or dataset-relative") from exc
        lines.append(f"sequence_reference = {reference_text}")
    if resolved.mirror:
        lines.append("mirror = true")
    if resolved.allow_unreviewed_backbone:
        lines.append("allow_unreviewed_backbone = true")
    # Always frozen: reopening a snapshot must not adopt a later default.
    lines.append(f"phaser_data = {resolved.phaser_data}")
    # Always freeze even an empty list: reopening a snapshot must not adopt a
    # later recipe version's phosphate sites.
    lines.append("allow_op3_sites = " + ", ".join(resolved.allow_op3_sites))
    if resolved.sequences:
        lines.extend(["", "[sequences]"])
        lines.extend(f"{chain} = {sequence}" for chain, sequence in resolved.sequences.items())
    if resolved.mutations:
        lines.extend(["", "[mutations]"])
        lines.extend(f"{site} = {ligand.token}" for site, ligand in resolved.mutations.items())
    if resolved.backbone_sites:
        lines.extend(["", "[backbones]"])
        lines.extend(f"{site} = {mode}" for site, mode in resolved.backbone_sites.items())
    return "\n".join(lines) + "\n"
