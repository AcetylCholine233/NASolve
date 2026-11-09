"""Immutable, portable campaign plans; no stage execution or workspace state.

DISCOVERED means only that input discovery, configuration, and model selection
succeeded. Scientific model, symmetry, chemistry, and external-tool preflights
remain the responsibility of the existing stage engines.
"""

from __future__ import annotations

import configparser
import copy
import hashlib
import json
import os
import re
import stat
import tempfile
from dataclasses import asdict, replace
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any

from .phosphate import PhosphateError, validate_op3_sites, validate_phosphate_intent
from .automr_input import (
    AutoMRInputError, AutoMRIntent, _parser, discover_dataset, locate_frames_directory,
    normalize_frame, read_intent, resolve_automr_input,
)
from .presets import PresetError, ProjectPreset, load_preset
from .campaign_threads import (
    SequenceThreadError, membership_by_dataset, parse_sequence_threads,
)


class CampaignError(RuntimeError):
    """A campaign cannot be created or its frozen state cannot be read safely."""


_STATE_DIRECTORY = "NASolveCampaign"
_SCOPE = "input_and_model_selection"
_EXCLUDED = {
    "automr", "postmr", "autosol", "autorefine", "cootgui", "refinedoctor",
    "nasolvecampaign", "env", "venv", "virtualenv", "conda", "miniconda",
    "miniconda3", "anaconda", "anaconda3", "__pycache__", "node_modules",
}
_HASH = re.compile(r"[0-9a-f]{64}\Z")
_AUTOMR_FIELDS = (
    "mode", "frame", "pair", "model", "sequence_file", "mirror", "allow_p1_standard",
)


def _canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _safe_relative(value: Any, label: str) -> str:
    if type(value) is not str or not value or "\x00" in value or "\\" in value:
        raise CampaignError(f"Invalid {label}: expected a safe relative path")
    posix, windows = PurePosixPath(value), PureWindowsPath(value)
    if (posix.is_absolute() or windows.drive or windows.root
            or any(part in {".", ".."} for part in value.split("/"))
            or str(posix) != value):
        raise CampaignError(f"Invalid {label}: path must remain relative to the campaign")
    return value


def _dataset_name(value: Any) -> str:
    name = _safe_relative(value, "dataset name")
    if "/" in name or name.startswith(".") or name.casefold() in _EXCLUDED:
        raise CampaignError(f"Invalid direct-child dataset name: {name!r}")
    return name


def _root(path: Path) -> Path:
    try:
        root = Path(path).expanduser().resolve(strict=True)
        if not root.is_dir():
            raise CampaignError(f"Campaign root is not a directory: {root}")
        return root
    except (OSError, RuntimeError, ValueError) as exc:
        raise CampaignError(f"Cannot access campaign root {path}: {exc}") from exc


def _contained(root: Path, path: Path, label: str) -> Path:
    try:
        resolved = path.resolve()
        resolved.relative_to(root)
    except (OSError, RuntimeError, ValueError) as exc:
        raise CampaignError(f"{label} escapes its declared anchor or cannot be resolved") from exc
    return resolved


def _file_identity(path: Path) -> tuple[str, int]:
    """Stream large reflection files and detect replacement during the read."""
    if not stat.S_ISREG(path.stat().st_mode):
        raise CampaignError(f"Expected a regular input file: {path.name}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise CampaignError(f"Expected a regular input file: {path.name}")
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
        after = os.fstat(handle.fileno())
    current = path.stat()
    identity = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns)
    if identity(before) != identity(after) or identity(after) != identity(current):
        raise CampaignError(f"Input changed while being read: {path.name}")
    return digest.hexdigest(), after.st_size


def _input_ref(root: Path, dataset: Path, path: Path) -> dict[str, Any]:
    relative = _safe_relative(path.relative_to(root).as_posix(), f"Input {path.name}")
    _contained(dataset, path, f"Input {path.name}")
    digest, size = _file_identity(path)
    return {
        "anchor": "campaign", "relative_path": relative,
        "sha256": digest, "size": size,
    }


def _resource_ref(staging: Path, data: bytes) -> dict[str, Any]:
    digest = _digest(data)
    relative = f"{_STATE_DIRECTORY}/resources/{digest}"
    target = staging / "resources" / digest
    if not target.exists():
        target.parent.mkdir(exist_ok=True)
        with target.open("xb") as handle:
            handle.write(data)
    return {"anchor": "campaign", "relative_path": relative, "sha256": digest, "size": len(data)}


def _resource_bytes(path: Path) -> bytes:
    digest, size = _file_identity(path)
    data = path.read_bytes()
    if len(data) != size or _digest(data) != digest:
        raise CampaignError(f"Resource changed while being read: {path.name}")
    return data

def _freeze_sequence_threads(
    root: Path, staging: Path,
) -> tuple[dict[str, Any] | None, dict[str, dict[str, object]]]:
    """Freeze optional root thread declarations and return explicit membership."""
    source = root / "nasolve-campaign.toml"
    if not source.exists() and not source.is_symlink():
        return None, {}
    if source.is_symlink() or not source.is_file():
        raise CampaignError("nasolve-campaign.toml must be a regular file in the campaign root")
    data = _resource_bytes(source)
    try:
        definitions = parse_sequence_threads(data)
    except SequenceThreadError as exc:
        raise CampaignError(str(exc)) from exc
    for thread in definitions.values():
        for dataset in thread["datasets"]:
            _dataset_name(dataset)
    frozen = {
        "source": _resource_ref(staging, data),
        "definitions": definitions,
    }
    return frozen, membership_by_dataset(definitions)



def _is_candidate_name(name: str) -> bool:
    path = Path(name)
    normalized = re.sub(r"[^a-z0-9]+", "", path.stem.casefold())
    return (path.suffix.casefold() == ".mtz"
            or (path.suffix.casefold() == ".cif" and normalized.startswith("data1"))
            or name.casefold() == "summary.html" or name == "nasolve.txt")


def _inventory(dataset: Path) -> list[str]:
    # Include dangling/escaping links as candidates so they receive a diagnostic.
    return sorted(path.name for path in dataset.iterdir()
                  if _is_candidate_name(path.name) and (path.is_file() or path.is_symlink()))


def _environment_directory(path: Path) -> bool:
    return (path / "pyvenv.cfg").is_file() or (path / "conda-meta").is_dir()


def _discover(root: Path, requested: tuple[str, ...] | None) -> list[Path]:
    if requested is not None:
        if not requested:
            raise CampaignError("Explicit dataset selection must not be empty")
        names = [_dataset_name(name) for name in requested]
        if len(set(names)) != len(names):
            raise CampaignError("Explicit dataset selection contains duplicate names")
        selected = []
        for name in sorted(names):
            path = root / name
            if path.is_symlink() or not path.is_dir() or _environment_directory(path):
                raise CampaignError(f"Requested dataset must be an existing direct-child directory: {name}")
            selected.append(path)
        return selected
    selected = []
    for path in sorted(root.iterdir()):
        if (path.name.startswith(".") or path.name.casefold() in _EXCLUDED
                or path.is_symlink() or not path.is_dir() or _environment_directory(path)):
            continue
        # Unreadable children are isolated as blocked records by _plan_dataset.
        try:
            candidate = bool(_inventory(path))
        except OSError:
            candidate = True
        if candidate:
            _dataset_name(path.name)
            selected.append(path)
    return selected


def _merged_intent(dataset: Path, preset: ProjectPreset) -> AutoMRIntent:
    config = dataset / "nasolve.txt"
    exists = config.exists() or config.is_symlink()
    intent = read_intent(config if exists else None)
    supplied: set[str] = set()
    if exists:
        parser = _parser()
        parser.read_string(config.read_text(encoding="utf-8"))
        supplied = set(parser["automr"])
    defaults = preset.automr_defaults
    for name in _AUTOMR_FIELDS:
        if name not in supplied and name in defaults:
            setattr(intent, name, defaults[name])
    return intent


def _intent_config(intent: AutoMRIntent) -> dict[str, Any]:
    return {
        "mode": intent.mode, "frame": intent.frame, "pair": intent.pair,
        "mirror": intent.mirror, "allow_p1_standard": intent.allow_p1_standard,
        "model_selector": intent.model,
        "model_family": intent.model_family,
        "sequence_reference": intent.sequence_reference,
        "sequences": dict(intent.sequences), "mutations": dict(intent.mutations),
        "allow_op3_sites": list(intent.allow_op3_sites),
        "backbones": dict(intent.backbone_sites),
        "allow_unreviewed_backbone": intent.allow_unreviewed_backbone,
    }


def _plan_dataset(root: Path, dataset: Path, preset: ProjectPreset, staging: Path,
                  frames_directory: Path | None,
                  sequence_thread: dict[str, object] | None = None) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "id": dataset.name, "relative_path": dataset.name,
        "status": "BLOCKED", "diagnostic": "", "effective_config": None,
        "inputs": {}, "discovery_inventory": [], "config_present": False,
        "duplicate_of": None, "duplicate_group": [],
    }
    try:
        inventory = _inventory(dataset)
        entry["discovery_inventory"] = inventory
        entry["config_present"] = "nasolve.txt" in inventory
        for name in inventory:
            entry["inputs"][f"discovery:{name}"] = _input_ref(root, dataset, dataset / name)
        files = discover_dataset(dataset)
        for role in ("reflections", "metadata", "summary"):
            path = getattr(files, role)
            entry["inputs"][role] = _input_ref(root, dataset, path)
        intent = _merged_intent(dataset, preset)
        dataset_sequence_reference = intent.sequence_reference
        if sequence_thread is not None and intent.sequence_reference is None:
            intent.sequence_reference = sequence_thread["sequence_reference"]
        entry["effective_config"] = _intent_config(intent)
        if intent.phaser_data not in {None, "F"}:
            raise AutoMRInputError(
                "phaser_data = I is not yet supported by campaigns (they use F); "
                "run 'nasolve automr' for this dataset or remove phaser_data"
            )
        if (not intent.mode or intent.mode.strip().casefold() != "standard"
                or not intent.frame or normalize_frame(intent.frame).name != "W"):
            raise AutoMRInputError("Campaign schema 1 supports only standard W/5W6W datasets")
        located_frames = locate_frames_directory(frames_directory, environ={})
        resolved = resolve_automr_input(dataset, intent, frames_dir=located_frames, environ={}, recipe=preset)
        applied_thread = None
        if sequence_thread is not None:
            applied_thread = {
                "id": sequence_thread["id"],
                "sequences": dict(sequence_thread["sequences"]),
                "site_codes": dict(sequence_thread["site_codes"]),
            }
            resolved = replace(resolved, sequence_thread=applied_thread)
        provider = resolved.model_provider
        if not isinstance(provider, dict):
            raise CampaignError("Resolved model provider provenance is missing")
        if provider.get("kind") == "standard-frame-catalogue":
            _contained(located_frames, resolved.model, "Selected catalogue model")
        elif provider.get("kind") == "explicit-standard-model":
            if provider.get("location") == "frame-catalogue":
                _contained(located_frames, resolved.model, "Explicit frame-catalogue model")
            elif provider.get("location") == "dataset":
                _contained(dataset, resolved.model, "Explicit dataset model")
            else:
                raise CampaignError("Explicit standard model has an unsupported provider location")
        else:
            raise CampaignError("Campaign schema 1 requires a standard-frame model provider")
        model_data = _resource_bytes(resolved.model)
        if not model_data:
            raise AutoMRInputError("Selected catalogue model is empty")
        model = _resource_ref(staging, model_data)
        entry["inputs"]["model"] = model
        sequence_reference = None
        if resolved.sequence_reference is not None:
            sequence_reference = _resource_ref(
                staging, _resource_bytes(resolved.sequence_reference)
            )
            entry["inputs"]["sequence_reference"] = sequence_reference
        sequence = resolved.frame_sequence_source
        if sequence is not None:
            _contained(located_frames, sequence, "Catalogue sequence")
            entry["inputs"]["frame_sequence"] = _resource_ref(staging, _resource_bytes(sequence))
        if intent.source is not None:
            entry["inputs"]["config"] = entry["inputs"]["discovery:nasolve.txt"]
        entry["effective_config"].update({
            "mode": resolved.mode, "frame": resolved.frame.name,
            "allow_op3_sites": list(resolved.allow_op3_sites),
            "phosphate_intent": resolved.phosphate_intent,
            "backbones": dict(resolved.backbone_sites),
            "allow_unreviewed_backbone": resolved.allow_unreviewed_backbone,
            "pair": resolved.pair_text, "pair_ligands": [asdict(item) for item in resolved.pair],
            "model": model, "model_name": resolved.model.name,
            "model_source": resolved.model_source,
            "model_selector": resolved.model_selector,
            "model_family": resolved.model_family,
            "model_provider": resolved.model_provider,
            "model_pair": (
                [asdict(item) for item in resolved.model_pair]
                if resolved.model_pair is not None else None
            ),
            "exact_pair_model": resolved.exact_pair_model,
            "catalogue_warnings": list(resolved.catalogue_warnings),
            "config_source": entry["inputs"].get("config"),
            "frame_sequence": entry["inputs"].get("frame_sequence"),
            "sequence_reference": resolved.sequence_reference_label,
            "sequence_reference_source": (
                "dataset" if dataset_sequence_reference is not None
                else "thread" if sequence_thread is not None else None
            ),
            "sequence_thread": applied_thread,
            "mutations": {site: asdict(ligand) for site, ligand in resolved.mutations.items()},
        })
        # Config parsing and discovery must describe the exact bytes frozen above.
        if inventory != _inventory(dataset):
            raise CampaignError("Discovery inputs changed during planning")
        for reference in entry["inputs"].values():
            if reference["relative_path"].startswith(dataset.name + "/"):
                digest, size = _file_identity(root / reference["relative_path"])
                if (digest, size) != (reference["sha256"], reference["size"]):
                    raise CampaignError("Dataset input changed during planning")
        entry["status"] = "DISCOVERED"
        entry["diagnostic"] = "Input and model selection complete; scientific stage preflights have not run."
    except (AutoMRInputError, CampaignError, OSError, UnicodeError, configparser.Error) as exc:
        entry["diagnostic"] = str(exc)
    return entry


def _freeze_preset(preset: ProjectPreset, staging: Path) -> dict[str, Any]:
    source = _resource_ref(staging, preset.source_bytes)
    if source["sha256"] != preset.sha256:
        raise CampaignError("Preset source snapshot does not match its checksum")
    policy = preset.to_dict()
    resources = {}
    for name, data in sorted(preset.resource_bytes.items()):
        resources[name] = _resource_ref(staging, data)
        expected = policy["resources"][name]
        if any(resources[name][key] != expected[key] for key in ("sha256", "size")):
            raise CampaignError(f"Preset resource snapshot does not match its checksum: {name}")
    for key in ("source", "sha256", "config_sha256"):
        policy.pop(key, None)
    return {
        "id": preset.id, "version": preset.version, "sha256": preset.sha256,
        "config_sha256": preset.config_sha256, "source": source,
        "policy": policy, "resources": resources,
    }


def _counts(entries: list[dict[str, Any]]) -> dict[str, int]:
    return {
        "total": len(entries),
        "discovered": sum(entry["status"] == "DISCOVERED" for entry in entries),
        "blocked": sum(entry["status"] == "BLOCKED" for entry in entries),
        "duplicate_groups": len({tuple(entry["duplicate_group"]) for entry in entries
                                 if entry["duplicate_group"]}),
    }


def _duplicates(entries: list[dict[str, Any]]) -> None:
    groups: dict[str, list[dict[str, Any]]] = {}
    for entry in entries:
        reflections = entry["inputs"].get("reflections")
        if reflections:
            groups.setdefault(reflections["sha256"], []).append(entry)
    for group in groups.values():
        if len(group) > 1:
            ids = [entry["id"] for entry in group]
            for index, entry in enumerate(group):
                entry["duplicate_group"] = ids
                entry["duplicate_of"] = ids[0] if index else None


def _publish(root: Path, staging: Path) -> Path:
    """Reserve exclusively, then publish the manifest with an atomic hard link.

    Portable rename can replace an existing empty directory. Exclusive mkdir
    plus hard-link publication refuses every existing destination and exposes
    plan.json only when all resources are present.
    """
    destination = root / _STATE_DIRECTORY
    owned: list[Path] = []
    reserved = False
    try:
        destination.mkdir()
        reserved = True
        resources = destination / "resources"
        resources.mkdir()
        owned.append(resources)
        for source in sorted((staging / "resources").iterdir()):
            target = resources / source.name
            os.link(source, target)
            owned.append(target)
        plan_path = destination / "plan.json"
        os.link(staging / "plan.json", plan_path)
        return plan_path
    except OSError:
        # Never recursively remove a destination that another process could use.
        for path in reversed(owned):
            try:
                path.rmdir() if path.is_dir() else path.unlink()
            except OSError:
                pass
        if reserved:
            try:
                destination.rmdir()
            except OSError:
                pass
        raise


def plan_campaign(root: Path, *, preset: str | Path = "5w6w",
                  datasets: tuple[str, ...] | None = None,
                  frames_directory: Path | None = None) -> dict[str, Any]:
    """Freeze direct-child datasets and reviewed resources without executing tools."""
    campaign = _root(root)
    destination = campaign / _STATE_DIRECTORY
    if destination.exists() or destination.is_symlink():
        raise CampaignError(f"Campaign state already exists; refusing to overwrite: {destination}")
    try:
        if frames_directory is not None:
            frames_directory = Path(frames_directory).expanduser().resolve(strict=True)
            if not frames_directory.is_dir():
                raise CampaignError("Explicit frames directory must be an existing directory")
        project = load_preset(preset)
        selected = _discover(campaign, datasets)
        if not selected:
            raise CampaignError("No dataset candidates found; choose a root containing dataset directories "
                                "or select an existing direct child explicitly")
        with tempfile.TemporaryDirectory(prefix=".nasolve-campaign-", dir=campaign) as temporary:
            staging = Path(temporary)
            frozen_preset = _freeze_preset(project, staging)
            frozen_threads, membership = _freeze_sequence_threads(campaign, staging)
            entries = [
                _plan_dataset(
                    campaign, path, project, staging, frames_directory,
                    membership.get(path.name),
                )
                for path in selected
            ]
            _duplicates(entries)
            payload = {
                "schema_version": 1, "state": "PLANNED", "validation_scope": _SCOPE,
                "preset": frozen_preset, "datasets": entries, "counts": _counts(entries),
                **({"sequence_threads": frozen_threads} if frozen_threads is not None else {}),
            }
            payload["fingerprint"] = _digest(_canonical(payload))
            _validate_plan(payload)
            with (staging / "plan.json").open("x", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True, indent=2, allow_nan=False)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            plan_path = _publish(campaign, staging)
        return {**payload, "plan_path": str(plan_path)}
    except (PresetError, SequenceThreadError, OSError, UnicodeError, RuntimeError, ValueError) as exc:
        if isinstance(exc, CampaignError):
            raise
        raise CampaignError(f"Could not plan campaign: {exc}") from exc


def _require(value: Any, kind: type, label: str) -> Any:
    if type(value) is not kind:
        raise CampaignError(f"Malformed campaign state: {label} must be {kind.__name__}")
    return value


def _required_fields(value: dict[str, Any], fields: set[str], label: str) -> None:
    missing = sorted(fields - value.keys())
    if missing:
        raise CampaignError(f"Malformed campaign state: {label} is missing field(s): {', '.join(missing)}")


def _validate_ref(value: Any, label: str, prefix: str | None = None) -> None:
    _require(value, dict, label)
    if value.get("anchor") != "campaign":
        raise CampaignError(f"Malformed campaign state: {label} must use the campaign anchor")
    relative = _safe_relative(value.get("relative_path"), label)
    if prefix is not None and not relative.startswith(prefix + "/"):
        raise CampaignError(f"Malformed campaign state: {label} is outside its declared input location")
    if type(value.get("sha256")) is not str or not _HASH.fullmatch(value["sha256"]):
        raise CampaignError(f"Malformed campaign state: {label} requires a SHA-256 checksum")
    if type(value.get("size")) is not int or value["size"] < 0:
        raise CampaignError(f"Malformed campaign state: {label} requires a nonnegative integer size")


def _validate_sequence_threads(
    payload: dict[str, Any], resource_prefix: str,
) -> dict[str, dict[str, object]]:
    value = payload.get("sequence_threads")
    if value is None:
        return {}
    value = _require(value, dict, "sequence_threads")
    if set(value) != {"source", "definitions"}:
        raise CampaignError("Malformed campaign state: invalid sequence_threads record")
    _validate_ref(value["source"], "sequence_threads.source", resource_prefix)
    definitions = _require(value["definitions"], dict, "sequence_threads.definitions")
    membership: dict[str, dict[str, object]] = {}
    for thread_id, thread in definitions.items():
        if (
            type(thread_id) is not str
            or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", thread_id) is None
        ):
            raise CampaignError("Malformed campaign state: invalid sequence thread id")
        thread = _require(thread, dict, f"sequence_threads.{thread_id}")
        if set(thread) != {"id", "datasets", "sequence_reference", "sequences", "site_codes"}:
            raise CampaignError("Malformed campaign state: invalid sequence thread fields")
        if thread.get("id") != thread_id or thread.get("sequence_reference") != "w-metal-scaffold":
            raise CampaignError("Malformed campaign state: inconsistent sequence thread identity")
        datasets = _require(thread.get("datasets"), list, f"{thread_id}.datasets")
        if datasets != sorted(set(datasets)) or not datasets:
            raise CampaignError("Malformed campaign state: invalid sequence thread membership")
        sequences = _require(thread.get("sequences"), dict, f"{thread_id}.sequences")
        site_codes = _require(thread.get("site_codes"), dict, f"{thread_id}.site_codes")
        if any(
            type(chain) is not str or len(chain) != 1
            or type(sequence) is not str or not sequence
            for chain, sequence in sequences.items()
        ):
            raise CampaignError("Malformed campaign state: invalid thread sequence overlay")
        if any(
            type(site) is not str or re.fullmatch(r"[A-Za-z0-9_]:[^:\\s]+", site) is None
            or type(code) is not str or re.fullmatch(r"[A-Z0-9]{1,5}", code) is None
            for site, code in site_codes.items()
        ):
            raise CampaignError("Malformed campaign state: invalid thread site chemistry")
        applied = {
            "id": thread_id,
            "sequences": sequences,
            "site_codes": site_codes,
        }
        for dataset in datasets:
            name = _dataset_name(dataset)
            if name in membership:
                raise CampaignError("Malformed campaign state: dataset belongs to multiple sequence threads")
            membership[name] = applied
    return membership


def _validate_plan(payload: Any) -> None:
    _require(payload, dict, "record")
    _required_fields(payload, {"schema_version", "state", "validation_scope", "preset",
                               "datasets", "counts", "fingerprint"}, "record")
    if type(payload.get("schema_version")) is not int or payload["schema_version"] != 1:
        raise CampaignError("Unsupported campaign schema_version; expected integer 1")
    if payload.get("state") != "PLANNED" or payload.get("validation_scope") != _SCOPE:
        raise CampaignError("Malformed campaign state: unsupported state or validation scope")
    preset = _require(payload.get("preset"), dict, "preset")
    _required_fields(preset, {"id", "version", "sha256", "config_sha256", "source",
                             "policy", "resources"}, "preset")
    for key in ("id", "version", "sha256", "config_sha256"):
        _require(preset.get(key), str, f"preset.{key}")
        if not preset[key] or (key.endswith("sha256") and not _HASH.fullmatch(preset[key])):
            raise CampaignError(f"Malformed campaign state: invalid preset.{key}")
    policy = _require(preset.get("policy"), dict, "preset.policy")
    _required_fields(policy, {"schema_version", "id", "version", "description", "automr",
                             "postmr", "autosol", "autorefine", "resources"}, "preset.policy")
    if (type(policy["schema_version"]) is not int or policy["schema_version"] != 1
            or policy["id"] != preset["id"] or policy["version"] != preset["version"]
            or _digest(_canonical(policy)) != preset["config_sha256"]):
        raise CampaignError("Malformed campaign state: inconsistent preset policy identity")
    _require(policy["description"], str, "preset.policy.description")
    for section in ("automr", "postmr", "autosol", "autorefine"):
        _require(policy[section], dict, f"preset.policy.{section}")
    if "chemistry" in policy:
        chemistry = _require(policy["chemistry"], dict, "preset.policy.chemistry")
        try:
            validate_op3_sites(chemistry.get("terminal_phosphate_sites", []))
        except PhosphateError as exc:
            raise CampaignError(f"Malformed campaign recipe chemistry: {exc}") from exc
    resource_prefix = f"{_STATE_DIRECTORY}/resources"
    _validate_ref(preset.get("source"), "preset.source", resource_prefix)
    if preset["source"]["sha256"] != preset["sha256"]:
        raise CampaignError("Malformed campaign state: inconsistent preset source checksum")
    resources = _require(preset.get("resources"), dict, "preset.resources")
    resource_policy = _require(policy["resources"], dict, "preset.policy.resources")
    if resources.keys() != resource_policy.keys():
        raise CampaignError("Malformed campaign state: inconsistent preset resource inventory")
    for name, reference in resources.items():
        _validate_ref(reference, f"preset.resources.{name}", resource_prefix)
        original = _require(resource_policy[name], dict, f"preset.policy.resources.{name}")
        _safe_relative(original.get("relative_path"), f"preset.policy.resources.{name}")
        if any(original.get(field) != reference[field] for field in ("sha256", "size")):
            raise CampaignError("Malformed campaign state: inconsistent preset resource identity")
    thread_membership = _validate_sequence_threads(payload, resource_prefix)
    entries = _require(payload.get("datasets"), list, "datasets")
    ids = []
    for entry in entries:
        _require(entry, dict, "dataset")
        _required_fields(entry, {"id", "relative_path", "status", "diagnostic", "effective_config",
                                  "inputs", "discovery_inventory", "config_present",
                                  "duplicate_of", "duplicate_group"}, "dataset")
        name = _dataset_name(entry.get("id"))
        if entry.get("relative_path") != name:
            raise CampaignError("Malformed campaign state: dataset id/path mismatch")
        ids.append(name)
        if entry.get("status") not in {"DISCOVERED", "BLOCKED"}:
            raise CampaignError("Malformed campaign state: invalid dataset status")
        _require(entry.get("diagnostic"), str, f"{name}.diagnostic")
        _require(entry.get("config_present"), bool, f"{name}.config_present")
        inventory = _require(entry.get("discovery_inventory"), list, f"{name}.discovery_inventory")
        for item in inventory:
            # Filenames here are an opaque discovery inventory, never operational
            # references. Unsafe filenames can therefore remain blocked records.
            if (type(item) is not str or "/" in item or "\x00" in item
                    or not _is_candidate_name(item)):
                raise CampaignError("Malformed campaign state: invalid discovery filename")
        if inventory != sorted(set(inventory)):
            raise CampaignError("Malformed campaign state: duplicate or unordered discovery inventory")
        if entry["config_present"] != ("nasolve.txt" in inventory):
            raise CampaignError("Malformed campaign state: inconsistent configuration presence")
        inputs = _require(entry.get("inputs"), dict, f"{name}.inputs")
        for role, reference in inputs.items():
            _validate_ref(reference, f"{name}.{role}",
                          resource_prefix if role in {"model", "frame_sequence", "sequence_reference"} else name)
        config = entry.get("effective_config")
        if config is not None:
            _require(config, dict, f"{name}.effective_config")
            for boolean in ("mirror", "allow_p1_standard"):
                _require(config.get(boolean), bool, f"{name}.{boolean}")
            if "allow_unreviewed_backbone" in config:
                _require(config.get("allow_unreviewed_backbone"), bool,
                         f"{name}.allow_unreviewed_backbone")
            if "backbones" in config:
                from .backbone import BackboneError, validate_backbone_sites
                try:
                    validate_backbone_sites(config["backbones"])
                except BackboneError as exc:
                    raise CampaignError(f"Malformed campaign backbone chemistry: {exc}") from exc
            for field in (
                "mode", "frame", "pair", "sequence_reference",
                "model_selector", "model_family",
            ):
                if config.get(field) is not None:
                    _require(config[field], str, f"{name}.{field}")
                    if field == "model_family" and re.fullmatch(
                        r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", config[field]
                    ) is None:
                        raise CampaignError(
                            f"Malformed campaign state: invalid {name}.model_family"
                        )
                    if field == "sequence_reference" and (
                        not config[field].strip()
                        or any(char in config[field] for char in "\r\n\x00")
                    ):
                        raise CampaignError(
                            f"Malformed campaign state: invalid {name}.sequence_reference"
                        )
            try:
                sites = validate_op3_sites(config.get("allow_op3_sites", []))
                if "phosphate_intent" in config:
                    intent = config["phosphate_intent"]
                    validate_phosphate_intent(intent, sites)
                    declaration = intent.get("recipe")
                    if declaration is not None:
                        for field in ("id", "version", "sha256", "config_sha256"):
                            if declaration[field] != preset[field]:
                                raise PhosphateError("Frozen chemistry and campaign recipe identity disagree")
                        expected = policy.get("chemistry", {}).get("terminal_phosphate_sites", [])
                        if declaration["terminal_phosphate_sites"] != expected or declaration["frame"] != config.get("frame"):
                            raise PhosphateError("Frozen chemistry and campaign recipe sites/frame disagree")
            except PhosphateError as exc:
                raise CampaignError(f"Malformed campaign OP3 request: {exc}") from exc
            for field in ("sequences", "mutations"):
                _require(config.get(field), dict, f"{name}.{field}")
            expected_thread = thread_membership.get(name)
            if config.get("sequence_thread") != expected_thread:
                raise CampaignError(
                    f"Malformed campaign state: inconsistent {name}.sequence_thread"
                )
            if "sequence_reference_source" in config:
                reference_source = config.get("sequence_reference_source")
                if reference_source not in {None, "dataset", "thread"}:
                    raise CampaignError(
                        f"Malformed campaign state: invalid {name}.sequence_reference_source"
                    )
                if reference_source == "thread":
                    if (
                        expected_thread is None
                        or config.get("sequence_reference")
                        != payload["sequence_threads"]["definitions"][expected_thread["id"]]["sequence_reference"]
                    ):
                        raise CampaignError(
                            f"Malformed campaign state: inconsistent {name} thread reference"
                        )
                elif reference_source == "dataset":
                    if config.get("sequence_reference") is None:
                        raise CampaignError(
                            f"Malformed campaign state: dataset reference source has no selector"
                        )
                elif config.get("sequence_reference") is not None:
                    raise CampaignError(
                        f"Malformed campaign state: selected reference has no provenance source"
                    )
            elif expected_thread is not None:
                raise CampaignError(
                    f"Malformed campaign state: thread-bound {name} lacks reference provenance"
                )
            for field, role in (("model", "model"), ("config_source", "config"),
                                ("frame_sequence", "frame_sequence")):
                if field in config and config[field] != inputs.get(role):
                    raise CampaignError(f"Malformed campaign state: inconsistent {name}.{field}")
            if entry["status"] == "DISCOVERED":
                selected_reference = config.get("sequence_reference")
                frozen_reference = inputs.get("sequence_reference")
                if (selected_reference is None) != (frozen_reference is None):
                    raise CampaignError(
                        f"Malformed campaign state: inconsistent {name}.sequence_reference"
                    )
        elif entry["status"] == "DISCOVERED":
            raise CampaignError("Malformed campaign state: discovered dataset requires effective configuration")
        if entry["status"] == "DISCOVERED":
            if any(role not in inputs for role in ("reflections", "metadata", "summary", "model")):
                raise CampaignError("Malformed campaign state: discovered dataset is missing frozen inputs")
            _required_fields(config, {"mode", "frame", "pair", "mirror", "allow_p1_standard",
                                      "model", "model_name", "model_source", "model_pair",
                                      "exact_pair_model", "catalogue_warnings", "config_source",
                                      "frame_sequence", "pair_ligands", "sequences", "mutations"},
                             f"{name}.effective_config")
            if "model_provider" in config:
                provider = _require(config["model_provider"], dict, f"{name}.model_provider")
                kind = provider.get("kind")
                if kind == "standard-frame-catalogue":
                    expected_provider_fields = {"kind", "selection", "frame", "location", "selector"}
                    selection = provider.get("selection")
                    selector = provider.get("selector")
                    expected_source = (
                        "standard frame catalogue (W; exact pair)"
                        if selection == "exact-pair"
                        else "standard frame catalogue (W; fallback C_G.pdb)"
                    )
                    if type(selector) is str:
                        selector_path = Path(_safe_relative(selector, f"{name}.model_provider.selector"))
                    else:
                        selector_path = Path("")
                    if (
                        set(provider) != expected_provider_fields
                        or selection not in {"exact-pair", "fallback"}
                        or provider.get("frame") != "W"
                        or provider.get("location") != "frame-catalogue"
                        or type(selector) is not str or not selector
                        or len(selector_path.parts) != 1
                        or selector_path.suffix.casefold() != ".pdb"
                        or config.get("model_selector") is not None
                        or config.get("model_family") is not None
                        or config.get("model_name") != selector
                        or config.get("model_pair") is None
                        or config.get("exact_pair_model") is not (selection == "exact-pair")
                        or config.get("model_source") != expected_source
                    ):
                        raise CampaignError(
                            f"Malformed campaign state: invalid {name}.model_provider"
                        )
                elif kind == "explicit-standard-model":
                    model_family = config.get("model_family")
                    expected_provider_fields = {
                        "kind", "selection", "frame", "location", "selector"
                    } | ({"construct_family"} if model_family is not None else set())
                    selector = provider.get("selector")
                    location = provider.get("location")
                    expected_source = (
                        f"explicit standard model (W; {location}:{selector})"
                        if isinstance(selector, str) and isinstance(location, str)
                        else None
                    )
                    if type(selector) is str:
                        selector_path = Path(_safe_relative(selector, f"{name}.model_provider.selector"))
                    else:
                        selector_path = Path("")
                    if (
                        set(provider) != expected_provider_fields
                        or provider.get("selection") != "user-forced"
                        or provider.get("frame") != "W"
                        or location not in {"dataset", "frame-catalogue"}
                        or type(selector) is not str or not selector
                        or selector_path.suffix.casefold() != ".pdb"
                        or selector_path.name != config.get("model_name")
                        or config.get("model_selector") != selector
                        or provider.get("construct_family") != model_family
                        or config.get("model_pair") is not None
                        or config.get("exact_pair_model") is not False
                        or config.get("model_source") != expected_source
                    ):
                        raise CampaignError(
                            f"Malformed campaign state: invalid {name}.model_provider"
                        )
                else:
                    raise CampaignError(
                        f"Malformed campaign state: unsupported {name}.model_provider"
                    )
            elif config.get("model_selector") is not None:
                raise CampaignError(
                    f"Malformed campaign state: model selector has no provider provenance"
                )
            if config["mode"] != "standard" or config["frame"] != "W" or not config["pair"]:
                raise CampaignError("Malformed campaign state: invalid discovered W configuration")
            _require(config["exact_pair_model"], bool, f"{name}.exact_pair_model")
            for field in ("model_name", "model_source"):
                _require(config[field], str, f"{name}.{field}")
            model_name_path = Path(_safe_relative(config["model_name"], f"{name}.model_name"))
            if len(model_name_path.parts) != 1 or model_name_path.suffix.casefold() != ".pdb":
                raise CampaignError(
                    f"Malformed campaign state: invalid {name}.model_name"
                )
            for field in ("pair_ligands", "catalogue_warnings"):
                _require(config[field], list, f"{name}.{field}")
            if config["model_pair"] is not None:
                _require(config["model_pair"], list, f"{name}.model_pair")
        group = _require(entry.get("duplicate_group"), list, f"{name}.duplicate_group")
        for member in group:
            _dataset_name(member)
        if entry.get("duplicate_of") is not None:
            _dataset_name(entry["duplicate_of"])
    if ids != sorted(set(ids)):
        raise CampaignError("Malformed campaign state: duplicate or unordered dataset ids")
    expected = copy.deepcopy(entries)
    for entry in expected:
        entry.update(duplicate_of=None, duplicate_group=[])
    _duplicates(expected)
    if any((left["duplicate_of"], left["duplicate_group"]) !=
           (right["duplicate_of"], right["duplicate_group"]) for left, right in zip(entries, expected)):
        raise CampaignError("Malformed campaign state: inconsistent duplicate groups")
    counts = _require(payload.get("counts"), dict, "counts")
    if any(type(value) is not int for value in counts.values()) or counts != _counts(entries):
        raise CampaignError("Malformed campaign state: inconsistent dataset counts")
    fingerprint = payload.get("fingerprint")
    unsigned = {key: value for key, value in payload.items() if key != "fingerprint"}
    if type(fingerprint) is not str or fingerprint != _digest(_canonical(unsigned)):
        raise CampaignError("Malformed campaign state: frozen inventory fingerprint mismatch")


def _integrity_issue(root: Path, reference: dict[str, Any]) -> str | None:
    relative = reference["relative_path"]
    try:
        path = _contained(root, root / relative, relative)
        boundary = (root / _STATE_DIRECTORY / "resources" if relative.startswith(_STATE_DIRECTORY + "/")
                    else root / PurePosixPath(relative).parts[0])
        _contained(boundary, path, relative)
        if not path.is_file():
            return f"Missing frozen file: {relative}"
        digest, size = _file_identity(path)
        if digest != reference["sha256"] or size != reference["size"]:
            return f"Changed frozen file: {relative} (SHA-256 or size mismatch)"
    except (CampaignError, OSError, RuntimeError, ValueError) as exc:
        return f"Cannot verify frozen file {relative}: {exc}"
    return None


def _json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise CampaignError(f"Malformed campaign state: duplicate JSON field {key!r}")
        result[key] = value
    return result


def campaign_status(root: Path) -> dict[str, Any]:
    """Read only the frozen collection, reporting drift without adoption or repair."""
    campaign = _root(root)
    plan_path = campaign / _STATE_DIRECTORY / "plan.json"
    try:
        _contained(campaign, plan_path, "Campaign state")
        with plan_path.open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_json_object)
        _validate_plan(payload)
    except (OSError, UnicodeError, ValueError, TypeError, RecursionError) as exc:
        raise CampaignError(f"Cannot read campaign plan: {exc}") from exc
    cache: dict[tuple[str, str, int], str | None] = {}

    def check(reference: dict[str, Any]) -> str | None:
        key = (reference["relative_path"], reference["sha256"], reference["size"])
        if key not in cache:
            cache[key] = _integrity_issue(campaign, reference)
        return cache[key]

    preset = payload["preset"]
    global_references = [preset["source"], *preset["resources"].values()]
    if payload.get("sequence_threads") is not None:
        global_references.append(payload["sequence_threads"]["source"])
    issues = [issue for reference in global_references if (issue := check(reference))]
    for entry in payload["datasets"]:
        local = [issue for reference in entry["inputs"].values() if (issue := check(reference))]
        dataset = campaign / entry["relative_path"]
        try:
            _contained(campaign, dataset, f"Dataset {entry['id']}")
            if dataset.is_symlink():
                local.append("Dataset directory became a symbolic link")
            else:
                current = _inventory(dataset)
                if current != entry["discovery_inventory"]:
                    added = sorted(set(current) - set(entry["discovery_inventory"]))
                    removed = sorted(set(entry["discovery_inventory"]) - set(current))
                    local.append(f"Discovery inventory changed; added={added}, removed={removed}")
        except (CampaignError, OSError, RuntimeError, ValueError) as exc:
            local.append(f"Cannot verify dataset directory: {exc}")
        entry["integrity_issues"] = sorted(set(local))
        entry["integrity"] = "DRIFT" if local else "OK"
    payload.update({
        "plan_path": str(plan_path), "integrity_issues": sorted(set(issues)),
        "integrity": "DRIFT" if issues or any(entry["integrity"] == "DRIFT"
                                               for entry in payload["datasets"]) else "OK",
    })
    return payload


__all__ = ["CampaignError", "plan_campaign", "campaign_status"]
