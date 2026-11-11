"""AutoFDP: measure anomalous-scatterer occupancy by f'' refinement.

PARTLY IMPLEMENTED. Step 1 works: ``prepare_start_model`` copies an AutoRefine
result, finds the metal HETATM sites and sets each one's B to the mean B of
the nucleotides around it (found with the NARestraints metal-pair search),
and records the copy as a new checkpoint. The f'' rounds and everything after
them are documented stubs that raise NotImplementedError. The full contract,
defaults and open questions are in docs/autofdp.md.

Method (supplements of the mmDNA papers, after Liu & Hendrickson 2013):
hold coordinates, B and occupancy; fix f' at its theoretical value for the
collection wavelength; refine f'' alone; adjust occupancy until refined f''
is within +/-1 e- of the theoretical f'' of the metal.
"""

from __future__ import annotations

import json
import math
import re
import shutil
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

from .checkpoints import (
    REUSABLE_STATUSES,
    CheckpointError,
    append_checkpoint,
    inherited_paths,
    next_checkpoint_id,
    read_registry,
    resolve_checkpoint,
)
from .model_assessment import file_sha256
from .run_context import artifact_reference

WORKFLOW = "autofdp"
SCHEMA_VERSION = 1

STATUS_PREPARED = "AUTOFDP_PREPARED"    # exit 0; step 1 only (B start)
STATUS_CONVERGED = "AUTOFDP_CONVERGED"  # exit 0
STATUS_REVIEW = "AUTOFDP_REVIEW"        # exit 2
STATUS_FAILED = "AUTOFDP_FAILED"        # exit 1

RECIPE_B_START = "AutoFDP/b-start"
RECIPE_FDP = "AutoFDP/fdp"
RECIPE_CONFIRM_B = "AutoFDP/confirm-b"
RECIPE_FINAL = "AutoFDP/final"

NOT_IMPLEMENTED_MESSAGE = (
    "this AutoFDP step is specified but not implemented yet; see docs/autofdp.md "
    "for what it will do."
)

#: Fallback radius (A) when a metal touches no base atom.
NEIGHBOUR_FALLBACK_DISTANCE = 6.0


class AutoFDPError(RuntimeError):
    """Raised when an AutoFDP run cannot be prepared or continued safely."""


@dataclass(frozen=True)
class AutoFDPSettings:
    """Command-line settings. Defaults marked (to confirm) in docs/autofdp.md."""

    sites: tuple[str, ...] = ()                       # empty: every anomalous scatterer
    elements: Mapping[str, str] = field(default_factory=dict)   # SITE -> ELEMENT
    candidates: tuple[str, ...] = ("AG", "HG")        # plus the model's elements
    theory: Mapping[str, float] = field(default_factory=dict)   # ELEMENT -> f''
    tolerance: float = 1.0                            # e-
    max_iterations: int = 8
    cycles: int = 3
    b_start: str = "surrounding"                      # surrounding | model | number
    b_step: float = 20.0                              # A^2
    occupancy_floor: float = 0.05
    map_sigma: float = 3.0
    confirm_b: bool = True
    final: bool = False
    plan_only: bool = False


@dataclass(frozen=True)
class SiteState:
    """One anomalous site in one AutoFDP round."""

    site: str                 # "B:101"
    atom_name: str
    element: str
    occupancy: float
    b_factor: float
    f_prime: float            # fixed at theory
    f_double_prime_theory: float
    f_double_prime_refined: float | None = None
    fofc_sigma: float | None = None
    anomalous_sigma: float | None = None
    decision: str | None = None


def plan_autofdp(
    run_directory: Path,
    mtz_dump_executable: Path,
    refine_executable: Path,
    settings: AutoFDPSettings,
    *,
    from_checkpoint: str | None = None,
    environment: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Preflight: check everything and return the plan, writing nothing.

    Checks that the run has a completed PostMR and a usable starting
    checkpoint; that the observations hold an anomalous quartet and a Free-R
    label (never regenerated); reads the wavelength from the MTZ header
    (AutoSol record as fallback, stop if neither); collects the sites (PostMR
    metal sites plus nucleotide heavy-atom candidates, or ``settings.sites``),
    each present exactly once in the starting model; computes theoretical
    f'/f'' for every site and candidate element (or uses ``settings.theory``);
    records the Phenix version; and checks with ``phenix.refine --dry_run``
    that the f''-only parameters are accepted. Raises AutoFDPError on any
    failure.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def execute_autofdp(
    run_directory: Path,
    refine_executable: Path,
    mtz_dump_executable: Path,
    *,
    phenix_version: str,
    settings: AutoFDPSettings,
    from_checkpoint: str | None = None,
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> dict[str, object]:
    """Run the whole calibration and return the AutoFDP report.

    1. ``plan_autofdp``.
    2. ``prepare_start_model`` (implemented): copy the AutoRefine result and
       set each metal's B from its neighbouring nucleotides, as a new
       checkpoint.
    3. Repeat up to ``settings.max_iterations``: ``fdp_round``, then
       ``decide_site_update`` for every unconverged site, and write the new
       occupancy/B as a model-edit checkpoint.
    4. ``identity_check`` after the first round.
    5. ``confirm_b_round`` unless disabled; ``final_round`` if requested.
    6. Write ``AutoFDP/fdp_NNN/report.json``, ``autofdp.log`` and
       ``fdp_table.tsv``. On convergence, bookmark the final checkpoint as
       ``autofdp-NNN``. Never change the current checkpoint.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def theoretical_scattering(
    elements: Sequence[str],
    wavelength: float,
    refine_executable: Path,
    environment: Mapping[str, str] | None = None,
    overrides: Mapping[str, float] | None = None,
) -> dict[str, dict[str, object]]:
    """Return {element: {f_prime, f_double_prime, source}} at the wavelength.

    Values come from the cctbx Henke table through ``phenix.python``, as in
    AutoRefine's calculated scattering; an override replaces f'' and is
    recorded with source ``user``. A non-finite or non-positive f'' is an
    AutoFDPError.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


@dataclass(frozen=True)
class Neighbour:
    """A nucleotide used for a metal's starting B."""

    residue: str          # "A:12"
    residue_name: str     # "DC"
    atoms: int            # atoms averaged (all atoms of the residue, first model)
    mean_b: float
    contact: str          # closest atom to the metal, e.g. "N4 2.45 A"


@dataclass(frozen=True)
class NeighbourB:
    """Mean B of a metal's neighbouring nucleotides and how they were chosen."""

    b_factor: float                 # atom-weighted mean over the neighbours
    method: str
    neighbours: tuple[Neighbour, ...]


@dataclass(frozen=True)
class BStartSite:
    """One metal site in AutoFDP step 1."""

    site: str             # "B:101"
    residue: str          # residue name, e.g. "AG"
    atom_name: str
    element: str
    alternate: str | None
    line_number: int      # 1-based line in the source model
    occupancy: float
    b_before: float
    b_after: float
    b_source: str         # "surrounding", "model" or "value"
    neighbour_b: NeighbourB | None


@dataclass(frozen=True)
class AutoFDPPreparation:
    """Result of AutoFDP step 1 (B start)."""

    run_directory: Path
    source_checkpoint: str
    source_model: Path                # the AutoRefine result that was copied
    sites: tuple[BStartSite, ...]
    warnings: tuple[str, ...]
    plan_only: bool
    output_directory: Path | None = None
    source_copy: Path | None = None
    start_model: Path | None = None
    checkpoint_id: str | None = None
    report_path: Path | None = None
    status: str = STATUS_PREPARED


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _metal_pairs():
    """The NARestraints metal-pair search (restraints.metal_pairs)."""
    try:
        from restraints import metal_pairs
    except ImportError as exc:
        raise AutoFDPError(
            "AutoFDP needs the NARestraints metal-pair search (restraints.metal_pairs), "
            "from the feature/anomalous-scatterers branch of NARestraints"
        ) from exc
    return metal_pairs


def _b_factor(line: str, line_number: int) -> float:
    try:
        value = float(line[60:66])
    except ValueError as exc:
        raise AutoFDPError(f"Line {line_number} has no readable B factor") from exc
    if not math.isfinite(value):
        raise AutoFDPError(f"Line {line_number} has a non-finite B factor")
    return value


def _occupancy(line: str, line_number: int) -> float:
    try:
        return float(line[54:60])
    except ValueError as exc:
        raise AutoFDPError(f"Line {line_number} has no readable occupancy") from exc


def surrounding_b(model_path: Path, metal_line: int) -> NeighbourB:
    """Mean B of the nucleotides around the metal atom on ``metal_line``.

    The neighbours are chosen with the NARestraints metal-pair search
    (restraints.metal_pairs), in this order:

    1. the two nucleotides of the metal-mediated base pair the ion bridges
       (``find_metal_mediated_pairs``: base N/O/S/Se within 3.6 A, two
       nucleotides that are not chain neighbours, closest pair kept);
    2. otherwise every nucleotide with a base N/O/S/Se within 3.6 A;
    3. otherwise every nucleotide with any atom within 6 A.

    The B factor is the mean over all atoms of those nucleotides (first model),
    rounded to 0.01. No nucleotide within 6 A is an AutoFDPError.
    """
    pairs_module = _metal_pairs()
    lines = model_path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms = pairs_module.read_atoms(model_path)
    metal = next((atom for atom in atoms if atom.line_number == metal_line), None)
    if metal is None:
        raise AutoFDPError(f"No atom on line {metal_line} of {model_path.name}")
    label = f"{metal.element} {pairs_module.residue_label(metal.residue)}"

    by_residue: dict[tuple[str, str], list[object]] = {}
    for atom in atoms:
        by_residue.setdefault(atom.residue, []).append(atom)
    nucleotides = {
        key for key, members in by_residue.items()
        if any(atom.name in ("C1'", "C1*") for atom in members)
    }

    def closest(predicate, limit: float) -> dict[tuple[str, str], tuple[float, str]]:
        found: dict[tuple[str, str], tuple[float, str]] = {}
        for atom in atoms:
            if atom.residue not in nucleotides or not predicate(atom):
                continue
            distance = math.dist(metal.xyz, atom.xyz)
            if distance <= limit and (
                atom.residue not in found or distance < found[atom.residue][0]
            ):
                found[atom.residue] = (distance, atom.name)
        return found

    pairs, _ = pairs_module.find_metal_mediated_pairs(model_path, metals=[metal.element])
    pair = next((item for item in pairs if item.metal.line_number == metal_line), None)
    if pair is not None:
        contacts = {
            contact.residue: (contact.distance, contact.atom_name)
            for contact in (pair.first, pair.second)
        }
        method = "bridged base pair (NARestraints metal-pair search)"
    else:
        contacts = closest(
            lambda atom: atom.name not in pairs_module.BACKBONE_ATOMS
            and atom.element in pairs_module.LIGAND_ELEMENTS,
            pairs_module.DEFAULT_MAX_DISTANCE,
        )
        method = f"base atoms within {pairs_module.DEFAULT_MAX_DISTANCE:g} A (no bridged pair)"
        if not contacts:
            contacts = closest(lambda atom: True, NEIGHBOUR_FALLBACK_DISTANCE)
            method = f"nucleotides within {NEIGHBOUR_FALLBACK_DISTANCE:g} A (no base contact)"
    if not contacts:
        raise AutoFDPError(
            f"No nucleotide within {NEIGHBOUR_FALLBACK_DISTANCE:g} A of {label}; remove "
            "the ion or leave it out with --site"
        )

    neighbours: list[Neighbour] = []
    every_b: list[float] = []
    for key, (distance, atom_name) in sorted(contacts.items(), key=lambda item: item[1][0]):
        members = by_residue[key]
        values = [_b_factor(lines[atom.line_number - 1], atom.line_number) for atom in members]
        every_b.extend(values)
        neighbours.append(Neighbour(
            residue=pairs_module.residue_label(key),
            residue_name=members[0].residue_name,
            atoms=len(values),
            mean_b=round(sum(values) / len(values), 2),
            contact=f"{atom_name} {distance:.2f} A",
        ))
    return NeighbourB(
        b_factor=round(sum(every_b) / len(every_b), 2),
        method=method,
        neighbours=tuple(neighbours),
    )


def parse_b_start(value: str | float) -> str | float:
    """``surrounding``, ``model`` or a B factor between 0 and 999.99."""
    if isinstance(value, str) and value.strip().casefold() in {"surrounding", "model"}:
        return value.strip().casefold()
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AutoFDPError(
            f"--b-start must be 'surrounding', 'model' or a number, not {value!r}"
        ) from exc
    if not math.isfinite(number) or not 0 < number < 1000:
        raise AutoFDPError("--b-start must be a B factor above 0 and below 1000")
    return number


def _metal_sites(model: Path, requested: Sequence[str]) -> tuple[list[tuple[dict, object]], list[str]]:
    from .postmr import PostMRPreparationError, scan_metal_anomalous_sites

    try:
        candidates, objects, warnings = scan_metal_anomalous_sites(model)
    except PostMRPreparationError as exc:
        raise AutoFDPError(str(exc)) from exc
    found = list(zip(candidates, objects))
    if not found:
        raise AutoFDPError(
            f"No metal HETATM (Z >= 21, outside nucleotides) in {model.name}"
        )
    if requested:
        wanted = [site.strip() for site in requested]
        known = {str(candidate["site"]) for candidate, _ in found}
        missing = [site for site in wanted if site not in known]
        if missing:
            raise AutoFDPError(
                f"No metal site {', '.join(missing)} in {model.name}; metal sites: "
                + ", ".join(sorted(known))
            )
        found = [item for item in found if str(item[0]["site"]) in wanted]
    return found, list(warnings)


def _set_b_factor(line: str, value: float) -> str:
    """Replace the B-factor columns (61-66) and keep the rest of the line."""
    if not 0 <= value < 1000:
        raise AutoFDPError(f"B factor {value} does not fit the PDB B column")
    body = line.rstrip("\r\n")
    ending = line[len(body):]
    body = body.ljust(66)
    return body[:60] + f"{value:6.2f}" + body[66:] + ending


def _next_fdp_directory(run: Path) -> Path:
    root = run / "AutoFDP"
    root.mkdir(parents=True, exist_ok=True)
    numbers = [
        int(match.group(1))
        for path in root.glob("fdp_*")
        if path.is_dir() and (match := re.fullmatch(r"fdp_(\d+)", path.name))
    ]
    destination = root / f"fdp_{max(numbers, default=0) + 1:03d}"
    destination.mkdir()
    return destination


def _file_reference(path: Path, run: Path) -> dict[str, object]:
    return {
        **artifact_reference(path, run),
        "sha256": file_sha256(path),
        "size": path.stat().st_size,
    }


def _site_record(site: BStartSite) -> dict[str, object]:
    record = {
        key: value for key, value in asdict(site).items()
        if key not in {"neighbour_b", "line_number"}
    }
    record["neighbour_method"] = site.neighbour_b.method if site.neighbour_b else None
    record["neighbour_mean_b"] = site.neighbour_b.b_factor if site.neighbour_b else None
    record["neighbours"] = (
        [asdict(item) for item in site.neighbour_b.neighbours] if site.neighbour_b else []
    )
    return record


def _log_lines(preparation: AutoFDPPreparation) -> list[str]:
    lines = [
        f"AutoFDP step 1 (B start), {_now()}",
        f"Source checkpoint: {preparation.source_checkpoint}",
        f"Source model: {preparation.source_model}",
    ]
    for site in preparation.sites:
        lines.append(
            f"{site.site} {site.atom_name} ({site.element}): occupancy {site.occupancy:.2f}; "
            f"B {site.b_before:.2f} -> {site.b_after:.2f} ({site.b_source})"
        )
        if site.neighbour_b is not None:
            lines.append(f"  neighbours: {site.neighbour_b.method}")
            for item in site.neighbour_b.neighbours:
                lines.append(
                    f"    {item.residue} {item.residue_name}: mean B {item.mean_b:.2f} "
                    f"over {item.atoms} atoms; closest {item.contact}"
                )
    lines.extend(f"Warning: {warning}" for warning in preparation.warnings)
    return lines


def prepare_start_model(
    run_directory: Path,
    *,
    from_checkpoint: str | None = None,
    sites: Sequence[str] = (),
    b_start: str | float = "surrounding",
    plan_only: bool = False,
) -> AutoFDPPreparation:
    """AutoFDP step 1: copy an AutoRefine result and set the metals' B.

    The starting checkpoint (``from_checkpoint``, default the current one)
    must be an AutoRefine refinement. Its model is copied unchanged to
    ``AutoFDP/fdp_NNN/source_model.pdb``; the metal HETATM sites are found as
    in PostMR (NARestraints ``find_heavy_atoms``: Z >= 21, outside nucleotides,
    not Li/Be/Na/Mg/Al/K/Ca), optionally restricted to ``sites``; and each
    site's B is set according to ``b_start``:

    - ``"surrounding"``: :func:`surrounding_b`, the mean B of the neighbouring
      nucleotides;
    - ``"model"``: unchanged;
    - a number: that B factor.

    Only the B column of those atom lines changes, in
    ``AutoFDP/fdp_NNN/start_model.pdb``, which is registered as a new
    model-edit checkpoint ``autofdp-NNN`` (kind ``manual``, recipe
    ``AutoFDP/b-start``, status REVIEW) under the source checkpoint. The
    current checkpoint is not changed. With ``plan_only`` nothing is written.
    """
    b_mode = parse_b_start(b_start)
    try:
        run, registry = read_registry(run_directory)
        item = resolve_checkpoint(registry, from_checkpoint)
    except CheckpointError as exc:
        raise AutoFDPError(str(exc)) from exc
    checkpoint = str(item.get("id"))
    if item.get("kind") != "refinement":
        refinements = [
            str(other.get("id")) for other in registry.get("checkpoints", [])
            if isinstance(other, Mapping) and other.get("kind") == "refinement"
        ]
        hint = (
            f"choose one with --from ({', '.join(refinements)})" if refinements
            else "run 'nasolve autorefine' first"
        )
        raise AutoFDPError(
            f"Checkpoint {checkpoint} is a {item.get('kind')} checkpoint, not an "
            f"AutoRefine result; {hint}"
        )
    if item.get("usable") is not True or item.get("status") not in REUSABLE_STATUSES:
        raise AutoFDPError(f"Checkpoint {checkpoint} is not reusable")
    try:
        model = inherited_paths(item, run)["model"]
    except CheckpointError as exc:
        raise AutoFDPError(str(exc)) from exc
    assert isinstance(model, Path)
    if model.suffix.lower() != ".pdb":
        raise AutoFDPError(f"AutoFDP needs a PDB model, not {model.name}")

    found, warnings = _metal_sites(model, sites)
    text = model.read_text(encoding="utf-8", errors="strict")
    lines = text.splitlines(keepends=True)
    prepared: list[BStartSite] = []
    for candidate, site_object in found:
        line_number = int(getattr(site_object, "line_number"))
        line = lines[line_number - 1] if 0 < line_number <= len(lines) else ""
        if not line.startswith("HETATM") or line[12:16].strip() != candidate["atom_name"]:
            raise AutoFDPError(
                f"Could not locate the {candidate['site']} {candidate['atom_name']} record "
                f"on line {line_number} of {model.name}"
            )
        if line_number < len(lines) and lines[line_number].startswith("ANISOU") and (
            lines[line_number][6:27] == line[6:27]
        ):
            raise AutoFDPError(
                f"{candidate['site']} {candidate['atom_name']} has anisotropic B (ANISOU); "
                "AutoFDP sets isotropic B only"
            )
        b_before = _b_factor(line, line_number)
        neighbour_b = None
        if b_mode == "surrounding":
            neighbour_b = surrounding_b(model, line_number)
            b_after, source = neighbour_b.b_factor, "surrounding"
        elif b_mode == "model":
            b_after, source = b_before, "model"
        else:
            b_after, source = round(float(b_mode), 2), "value"
        prepared.append(BStartSite(
            site=str(candidate["site"]),
            residue=str(candidate.get("residue") or ""),
            atom_name=str(candidate["atom_name"]),
            element=str(candidate["element"]),
            alternate=candidate.get("alternate") or None,
            line_number=line_number,
            occupancy=_occupancy(line, line_number),
            b_before=b_before,
            b_after=b_after,
            b_source=source,
            neighbour_b=neighbour_b,
        ))

    preparation = AutoFDPPreparation(
        run_directory=run,
        source_checkpoint=checkpoint,
        source_model=model,
        sites=tuple(prepared),
        warnings=tuple(warnings),
        plan_only=plan_only,
    )
    if plan_only:
        return preparation

    destination = _next_fdp_directory(run)
    source_copy = destination / "source_model.pdb"
    shutil.copyfile(model, source_copy)
    if file_sha256(source_copy) != file_sha256(model):
        raise AutoFDPError("The copy of the AutoRefine model does not match its source")
    edited = list(lines)
    for site in prepared:
        edited[site.line_number - 1] = _set_b_factor(lines[site.line_number - 1], site.b_after)
    start_model = destination / "start_model.pdb"
    start_model.write_text("".join(edited), encoding="utf-8", newline="")

    report_path = destination / "report.json"
    checkpoint_id = next_checkpoint_id(registry, "autofdp")
    edits = [
        {"site": site.site, "atom_name": site.atom_name, "b_before": site.b_before,
         "b_after": site.b_after, "b_source": site.b_source}
        for site in prepared
    ]
    payload = {
        "id": checkpoint_id,
        "parent": checkpoint,
        # A model edit outside phenix.refine, like an imported model, so the
        # viewer shows it with its parent's maps and AutoRefine accepts it.
        "kind": "manual",
        "status": "REVIEW",
        "usable": True,
        "recipe": RECIPE_B_START,
        "label": None,
        "model": _file_reference(start_model, run),
        "observations": item.get("observations"),
        "phases": item.get("phases"),
        "restraints": item.get("restraints", []),
        **(
            {"terminal_geometry_protection": item["terminal_geometry_protection"]}
            if "terminal_geometry_protection" in item else {}
        ),
        "metrics": {},
        "compatibility": {
            "validated": False,
            "source": f"AutoFDP B start from {checkpoint}",
            "review_required": True,
        },
        "autofdp": {
            "step": "b-start",
            "record": artifact_reference(destination, run),
            "edits": edits,
        },
    }
    preparation = replace(
        preparation,
        output_directory=destination,
        source_copy=source_copy,
        start_model=start_model,
        checkpoint_id=checkpoint_id,
        report_path=report_path,
    )
    report: dict[str, object] = {
        "workflow": WORKFLOW,
        "schema_version": SCHEMA_VERSION,
        "status": STATUS_PREPARED,
        "created_utc": _now(),
        "steps": {"b_start": "done", "fdp_rounds": "not implemented"},
        "source_checkpoint": checkpoint,
        "source_model": {
            **_file_reference(source_copy, run),
            "copied_from": artifact_reference(model, run),
        },
        "start_model": _file_reference(start_model, run),
        "start_checkpoint": checkpoint_id,
        "settings": {"b_start": b_mode, "sites": list(sites)},
        "sites": [_site_record(site) for site in prepared],
        "warnings": list(warnings),
    }
    try:
        append_checkpoint(run, registry, payload, select=False)
    except CheckpointError as exc:
        report.update(status=STATUS_FAILED, start_checkpoint=None, error=str(exc))
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        raise AutoFDPError(f"Could not register the AutoFDP checkpoint: {exc}") from exc
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    (destination / "autofdp.log").write_text(
        "\n".join(_log_lines(preparation)) + "\n", encoding="utf-8"
    )
    return preparation


def fdp_round(
    run_directory: Path,
    checkpoint: str,
    sites: Sequence[SiteState],
    settings: AutoFDPSettings,
    refine_executable: Path,
    mtz_dump_executable: Path,
    *,
    phenix_version: str,
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> tuple[str, list[SiteState]]:
    """One f''-only refinement, recorded as an AutoRefine round.

    Strategy ``group_anomalous`` only (no coordinates, B or occupancy for any
    atom); one anomalous group per site with f' fixed at theory and f''
    refined from 0; ``settings.cycles`` macrocycles; recipe ``AutoFDP/fdp``.
    Returns the new checkpoint ID and the sites with refined f'' and the
    mFo-DFc and anomalous map values at each site (``map_value_at_site``).
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def map_value_at_site(
    map_coefficients: Path,
    labels: str,
    coordinates: tuple[float, float, float],
    phenix_python: Path,
    environment: Mapping[str, str] | None = None,
) -> float:
    """Value of a map (mFo-DFc or anomalous) at a point, in sigma units."""
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def decide_site_update(site: SiteState, previous: SiteState | None,
                       settings: AutoFDPSettings) -> SiteState:
    """Apply the docs/autofdp.md step-3 rules to one site.

    Converged if |refined - theory| <= tolerance. Otherwise proposed
    occupancy q' = q * refined / theory; adjust B by ``b_step`` instead when
    the mFo-DFc value at the site contradicts the occupancy change; cap at
    1.0 (flag "f'' above theory at full occupancy"); below
    ``occupancy_floor`` send to identity review; flag "unresponsive" when f''
    barely changed after an occupancy change. Returns the site with its next
    occupancy, B and decision.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def identity_check(
    site: SiteState,
    theory: Mapping[str, Mapping[str, object]],
    settings: AutoFDPSettings,
) -> dict[str, object]:
    """Which candidate elements explain the refined f''.

    For each candidate E, implied occupancy = q * refined / theory(E); E is
    plausible if that is above the floor and at most 1.0. Reports
    consistent, ambiguous (several plausible) or suggest-other (the model's
    element is not plausible). Never changes the element.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def confirm_b_round(
    run_directory: Path,
    checkpoint: str,
    sites: Sequence[SiteState],
    settings: AutoFDPSettings,
    refine_executable: Path,
    mtz_dump_executable: Path,
    *,
    phenix_version: str,
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> tuple[str, list[SiteState]]:
    """Refine only the sites' B with occupancy, f' and f'' fixed.

    Recipe ``AutoFDP/confirm-b``. A B change above 25% sends the site to
    review.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def final_round(
    run_directory: Path,
    checkpoint: str,
    sites: Sequence[SiteState],
    settings: AutoFDPSettings,
    refine_executable: Path,
    mtz_dump_executable: Path,
    *,
    phenix_version: str,
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> tuple[str, list[SiteState]]:
    """Optional Adv. Mater. ending: XYZ, B, occupancy and f'' together.

    f' stays fixed at theory. Recipe ``AutoFDP/final``. Reports how far
    occupancy and f'' move.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


__all__ = [
    "AutoFDPError",
    "AutoFDPPreparation",
    "AutoFDPSettings",
    "BStartSite",
    "Neighbour",
    "NeighbourB",
    "SiteState",
    "NOT_IMPLEMENTED_MESSAGE",
    "STATUS_CONVERGED",
    "STATUS_PREPARED",
    "STATUS_FAILED",
    "STATUS_REVIEW",
    "WORKFLOW",
    "confirm_b_round",
    "decide_site_update",
    "execute_autofdp",
    "fdp_round",
    "final_round",
    "identity_check",
    "map_value_at_site",
    "parse_b_start",
    "plan_autofdp",
    "prepare_start_model",
    "surrounding_b",
    "theoretical_scattering",
]
