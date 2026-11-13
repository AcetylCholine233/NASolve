"""AutoFDP: measure metal-site occupancy by f'' refinement.

PARTLY IMPLEMENTED. AutoFDP needs only a starting model (PDB) and an MTZ file;
it reads no NASolve run report or checkpoint.

- From the MTZ (``phenix.mtz.dump``): the anomalous observations
  (F(+)/SIGF(+)/F(-)/SIGF(-), or the intensity equivalent), the Free-R flags
  (never regenerated) and the wavelength (MTZ header; ``--wavelength`` if the
  header has none).
- From the model: the metal HETATM sites (NARestraints heavy-atom search, as
  in PostMR) and the B factors of the nucleotides around each one.

Implemented:

- Step 1 (``prepare_start_model``): copy the model and set each metal's B to
  the mean B of the nucleotides around it (NARestraints metal-pair search).
- Step 2: ``run_readyset`` (phenix.ready_set on that copy, no hydrogens, for
  the restraint dictionaries of non-standard residues; curated NASolve
  dictionaries and ``--cif`` files are given to it), then ``fdp_round``: one
  phenix.refine round on ReadySet's model with every strategy off except
  ``group_anomalous`` (f' fixed at theory, f'' refined from 0), and each
  refined f'' compared with its theoretical value.
- Step 3 (occupancy adjustment): the site furthest from theory first, one
  site at a time, its occupancy changed (``propose_occupancy``: lower when
  f'' is below theory, higher when above) and the same f''-only round run
  again, until its f'' is within the tolerance or it has had
  ``max_iterations`` adjustments; then the next site.

``plan_autofdp`` is the preflight; ``execute_autofdp`` runs steps 1 to 3.
The identity check, B confirmation and final round are documented stubs that
raise NotImplementedError. The full contract, defaults and open questions are
in docs/autofdp.md.

Method (supplements of the mmDNA papers, after Liu & Hendrickson 2013):
hold coordinates, B and occupancy; fix f' at its theoretical value for the
collection wavelength; refine f'' alone; adjust occupancy until refined f''
is within +/-1 e- of the theoretical f'' of the metal.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Mapping, Sequence

from .model_assessment import file_sha256
from .run_context import artifact_reference

WORKFLOW = "autofdp"
# 2: standalone (a model and an MTZ). 1 was the AutoRefine-checkpoint layout.
SCHEMA_VERSION = 2

STATUS_PREPARED = "AUTOFDP_PREPARED"           # step 1 written, round not yet run
STATUS_ROUND_DONE = "AUTOFDP_FDP_ROUND_DONE"    # one f''-only round finished
STATUS_CONVERGED = "AUTOFDP_CONVERGED"          # exit 0: every site within tolerance
STATUS_REVIEW = "AUTOFDP_REVIEW"                # exit 2: a site did not converge
STATUS_FAILED = "AUTOFDP_FAILED"                # exit 1

#: One occupancy adjustment changes a site's occupancy at most this many fold.
OCCUPANCY_MAX_RATIO = 4.0
#: Default number of occupancy adjustments per site (``--max-iterations``).
DEFAULT_MAX_ITERATIONS = 5

NOT_IMPLEMENTED_MESSAGE = (
    "this AutoFDP step is specified but not implemented yet; see docs/autofdp.md "
    "for what it will do."
)

#: Fallback radius (A) when a metal touches no base atom.
NEIGHBOUR_FALLBACK_DISTANCE = 6.0

MODEL_SUFFIXES = (".pdb", ".ent")
ANOMALOUS_LABEL_SETS = (
    ("F(+)", "SIGF(+)", "F(-)", "SIGF(-)"),
    ("I(+)", "SIGI(+)", "I(-)", "SIGI(-)"),
)
#: Plausible X-ray wavelengths (A), as for the MTZ-header reader.
WAVELENGTH_RANGE = (0.1, 5.0)


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
    max_iterations: int = DEFAULT_MAX_ITERATIONS     # occupancy adjustments per site
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


def parse_theory_overrides(values: Sequence[str]) -> dict[str, float]:
    """``ELEMENT=FDP`` strings (from ``--theory``) as {ELEMENT: f''}."""
    overrides: dict[str, float] = {}
    for text in values:
        element, separator, number = str(text).partition("=")
        element = element.strip().upper()
        if not separator or re.fullmatch(r"[A-Z]{1,2}", element) is None:
            raise AutoFDPError(f"--theory must look like HG=10.17, not {text!r}")
        try:
            value = float(number)
        except ValueError as exc:
            raise AutoFDPError(f"--theory {text!r}: f'' must be a number") from exc
        if not math.isfinite(value) or value <= 0:
            raise AutoFDPError(f"--theory {text!r}: f'' must be positive")
        if element in overrides and overrides[element] != value:
            raise AutoFDPError(f"--theory gives two values for {element}")
        overrides[element] = value
    return overrides


def theoretical_scattering(
    elements: Sequence[str],
    wavelength: float,
    refine_executable: Path,
    environment: Mapping[str, str] | None = None,
    overrides: Mapping[str, float] | None = None,
) -> dict[str, dict[str, object]]:
    """Return {element: {f_prime, f_double_prime, f_double_prime_henke, source}}.

    Values come from the cctbx Henke table at ``wavelength`` through the
    ``phenix.python`` beside ``refine_executable`` (as in AutoRefine). An
    override (``--theory``) replaces an element's f'' and is recorded with
    source ``user (--theory)``; f' always stays the Henke value.
    """
    wanted = sorted({element.upper() for element in elements})
    if not wanted:
        raise AutoFDPError("No element to calculate scattering for")
    overrides = {key.upper(): value for key, value in (overrides or {}).items()}
    unused = sorted(set(overrides) - set(wanted))
    if unused:
        raise AutoFDPError(
            f"--theory names {', '.join(unused)}, which is not an anomalous scatterer "
            f"here ({', '.join(wanted)})"
        )
    script = (
        "import json\nfrom cctbx.eltbx import henke\nvalues = {}\n"
        f"for element in {wanted!r}:\n"
        f"    value = henke.table(element).at_angstrom({float(wavelength)!r})\n"
        "    values[element] = [value.fp(), value.fdp()]\n"
        "print('NASOLVE_SCATTERING_JSON:' + json.dumps(values))\n"
    )
    python = refine_executable.expanduser().resolve().parent / "phenix.python"
    try:
        completed = subprocess.run(
            [str(python), "-c", script],
            env=dict(environment) if environment is not None else None,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=False, timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise AutoFDPError(f"Cannot calculate anomalous scattering with {python}: {exc}") from exc
    matches = re.findall(r"^NASOLVE_SCATTERING_JSON:(\{.*\})$", completed.stdout, re.M)
    if completed.returncode or len(matches) != 1:
        raise AutoFDPError(
            f"phenix.python returned no scattering values for {', '.join(wanted)}"
        )
    values = json.loads(matches[0])
    result: dict[str, dict[str, object]] = {}
    for element in wanted:
        pair = values.get(element) if isinstance(values, dict) else None
        if not isinstance(pair, list) or len(pair) != 2 or any(
            isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
            for v in pair
        ):
            raise AutoFDPError(f"Invalid calculated scattering values for {element}")
        theory = overrides.get(element, float(pair[1]))
        if not math.isfinite(theory) or theory <= 0:
            raise AutoFDPError(
                f"The theoretical f'' of {element} at {wavelength:g} A is {theory:g}; "
                "give one with --theory"
            )
        result[element] = {
            "f_prime": float(pair[0]),
            "f_double_prime": theory,
            "f_double_prime_henke": float(pair[1]),
            "source": "user (--theory)" if element in overrides else "cctbx-henke",
        }
    return result


@dataclass(frozen=True)
class ReflectionInputs:
    """What AutoFDP takes from the MTZ file."""

    mtz: Path
    observation_labels: tuple[str, ...]
    free_r_label: str
    free_r_test_value: int
    wavelength: float
    wavelength_source: Mapping[str, object]
    resolution_limit: float | None


def read_reflections(
    mtz: Path,
    mtz_dump_executable: Path,
    environment: Mapping[str, str] | None = None,
    *,
    wavelength: float | None = None,
) -> ReflectionInputs:
    """Read labels, Free-R flags and wavelength from ``phenix.mtz.dump``.

    The anomalous array is F(+)/SIGF(+)/F(-)/SIGF(-) (or I(+)/SIGI(+)/I(-)/
    SIGI(-)); without one, f'' cannot be refined and AutoFDP stops. The Free-R
    label must exist (flags are never generated); its test value is 1 for
    Phenix ``R-free-flags`` and 0 for the CCP4 labels (FreeR_flag, FREE,
    FreeR), as AutoRefine uses for STARANISO data. The wavelength is the MTZ
    header value of the dataset holding the data, unless ``wavelength`` is
    given.
    """
    from .autorefine import (
        AutoRefineError,
        _first_label_set,
        _free_r_label,
        _resolution_limit,
        _run_mtz_dump,
    )
    from .reflection_data import ReflectionDataError, mtz_wavelength

    mtz = mtz.expanduser().resolve()
    if not mtz.is_file():
        raise AutoFDPError(f"No MTZ file: {mtz}")
    try:
        dump = _run_mtz_dump(mtz, mtz_dump_executable, environment)
        free_r = _free_r_label(dump)
    except AutoRefineError as exc:
        raise AutoFDPError(str(exc)) from exc
    labels = _first_label_set(dump, ANOMALOUS_LABEL_SETS)
    if labels is None:
        raise AutoFDPError(
            f"{mtz.name} has no anomalous data (F(+)/SIGF(+)/F(-)/SIGF(-) or "
            "I(+)/SIGI(+)/I(-)/SIGI(-)); f'' cannot be refined"
        )
    try:
        header, dataset = mtz_wavelength(dump)
        header_error = None
    except ReflectionDataError as exc:
        header, dataset, header_error = None, None, str(exc)
    if wavelength is not None:
        low, high = WAVELENGTH_RANGE
        if not math.isfinite(wavelength) or not low <= wavelength <= high:
            raise AutoFDPError(f"--wavelength must be between {low:g} and {high:g} A")
        value = float(wavelength)
        source = {"method": "--wavelength", "mtz_header": header, "mtz_error": header_error}
        if header is not None and abs(header - value) > 1e-4:
            source["warning"] = (
                f"--wavelength {value:g} A differs from the MTZ header ({header:g} A)"
            )
    elif header is not None:
        value = header
        source = {"method": "MTZ header (phenix.mtz.dump)", "mtz_dataset": dataset}
    else:
        raise AutoFDPError(f"{header_error}; give the wavelength with --wavelength")
    return ReflectionInputs(
        mtz=mtz,
        observation_labels=tuple(labels),
        free_r_label=free_r,
        free_r_test_value=1 if free_r == "R-free-flags" else 0,
        wavelength=value,
        wavelength_source=source,
        resolution_limit=_resolution_limit(dump),
    )


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
class BStart:
    """AutoFDP step 1 (B start): the metal sites and their B before and after."""

    model: Path                       # the starting model given on the command line
    sites: tuple[BStartSite, ...]
    warnings: tuple[str, ...]
    source_copy: Path | None = None   # unchanged copy, once written
    start_model: Path | None = None   # copy with the metal B set, once written


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


def checked_model(path: Path) -> Path:
    """The starting model: an existing PDB-format file."""
    model = path.expanduser().resolve()
    if not model.is_file():
        raise AutoFDPError(f"No model file: {model}")
    if model.suffix.lower() not in MODEL_SUFFIXES:
        raise AutoFDPError(f"AutoFDP needs a PDB-format model (.pdb), not {model.name}")
    return model


def prepare_start_model(
    model: Path,
    *,
    sites: Sequence[str] = (),
    b_start: str | float = "surrounding",
    output_directory: Path | None = None,
) -> BStart:
    """AutoFDP step 1: find the metal sites and set their B.

    The metal HETATM sites are found as in PostMR (NARestraints
    ``find_heavy_atoms``: Z >= 21, outside nucleotides, not Li/Be/Na/Mg/Al/
    K/Ca), optionally restricted to ``sites``, and each site's B is set
    according to ``b_start``:

    - ``"surrounding"``: :func:`surrounding_b`, the mean B of the neighbouring
      nucleotides;
    - ``"model"``: unchanged;
    - a number: that B factor.

    Without ``output_directory`` nothing is written (the plan). With it, the
    model is copied unchanged to ``source_model.pdb`` and, with only the B
    column of the metal lines changed, to ``start_model.pdb``.
    """
    model = checked_model(model)
    b_mode = parse_b_start(b_start)
    found, warnings = _metal_sites(model, sites)
    try:
        text = model.read_text(encoding="utf-8", errors="strict")
    except (OSError, UnicodeDecodeError) as exc:
        raise AutoFDPError(f"Could not read {model.name}: {exc}") from exc
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
    result = BStart(model=model, sites=tuple(prepared), warnings=tuple(warnings))
    if output_directory is None:
        return result

    output_directory.mkdir(parents=True, exist_ok=True)
    source_copy = output_directory / "source_model.pdb"
    shutil.copyfile(model, source_copy)
    if file_sha256(source_copy) != file_sha256(model):
        raise AutoFDPError("The copy of the starting model does not match its source")
    edited = list(lines)
    for site in prepared:
        edited[site.line_number - 1] = _set_b_factor(lines[site.line_number - 1], site.b_after)
    start_model = output_directory / "start_model.pdb"
    start_model.write_text("".join(edited), encoding="utf-8", newline="")
    return replace(result, source_copy=source_copy, start_model=start_model)


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


def _site_selection(site: str, atom_name: str, alternate: str | None) -> str:
    chain, residue = site.split(":", 1)
    selection = f"chain {chain} and resid {residue} and name {atom_name}"
    if alternate:
        selection += f" and altloc {alternate}"
    return selection


def _site_key(site: str, atom_name: str, alternate: str | None) -> tuple[str, ...]:
    """The _coordinate_records key of a site such as "B:101" or "B:101A"."""
    chain, residue = site.split(":", 1)
    match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", residue)
    resseq, icode = (match.group(1), match.group(2)) if match else (residue, "")
    return (chain if chain != "_" else "", resseq, icode, atom_name, alternate or "")


def _record_key(line: str) -> tuple[str, ...]:
    """(chain, resseq, icode, name, altloc) of an ATOM/HETATM line."""
    return (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(),
            line[12:16].strip(), line[16:17].strip())


def _coordinate_records(path: Path) -> dict[tuple[str, ...], tuple[float, ...]]:
    """(chain, resseq, icode, name, altloc) -> (x, y, z, occupancy, B)."""
    records: dict[tuple[str, ...], tuple[float, ...]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = _record_key(line)
        try:
            records[key] = tuple(float(line[a:b]) for a, b in
                                 ((30, 38), (38, 46), (46, 54), (54, 60), (60, 66)))
        except ValueError:
            continue
    return records


def held_parameter_changes(start: Path, refined: Path) -> dict[str, object]:
    """Largest coordinate, occupancy and B change between two models.

    An f''-only round must leave all three alone; this is the check.
    """
    before, after = _coordinate_records(start), _coordinate_records(refined)
    common = before.keys() & after.keys()
    shift = max((math.dist(before[key][:3], after[key][:3]) for key in common), default=0.0)
    occupancy = max((abs(before[key][3] - after[key][3]) for key in common), default=0.0)
    b_factor = max((abs(before[key][4] - after[key][4]) for key in common), default=0.0)
    return {
        "atoms_compared": len(common),
        "atoms_only_in_start": len(before.keys() - after.keys()),
        "atoms_only_in_refined": len(after.keys() - before.keys()),
        "max_coordinate_shift": round(shift, 4),
        "max_occupancy_change": round(occupancy, 4),
        "max_b_change": round(b_factor, 3),
    }


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


@dataclass(frozen=True)
class ReadySetResult:
    """AutoFDP's ReadySet step: restraint dictionaries for the f'' round."""

    model: Path                       # the model phenix.refine gets (B from step 1 kept)
    dictionaries: tuple[Path, ...]    # restraint CIFs for phenix.refine
    generated: Path | None            # ReadySet's <model>.ligands.cif
    input_dictionaries: Path | None   # curated + --cif dictionaries given to ReadySet
    curated_codes: tuple[str, ...]    # NASolve-curated residues found in the model
    log: Path
    command: tuple[str, ...]
    mode: str                         # "updated model" or "no unknown residues"
    changes: Mapping[str, object]     # start model -> ReadySet model (must be none)


def _component_codes(path: Path) -> set[str]:
    """Residue codes a restraint CIF defines atoms for."""
    from .ligand_profiles import _blocks, _component_ids

    codes: set[str] = set()
    for block in _blocks(path):
        codes |= _component_ids(block)
    return codes


def _hydrogen_count(path: Path) -> int:
    count = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            element = line[76:78].strip().upper() if len(line) > 76 else ""
            name = line[12:16].strip().upper()
            if element in {"H", "D"} or (not element and name[:1] in {"H", "D"}
                                          and line[17:20].strip() != "HG"):
                count += 1
    return count


def run_readyset(
    start_model: Path,
    directory: Path,
    ready_set_executable: Path,
    *,
    extra_dictionaries: Sequence[Path] = (),
    environment: Mapping[str, str] | None = None,
) -> ReadySetResult:
    """Run phenix.ready_set on the step-1 model for the f'' round's restraints.

    The actions are PostMR's (``postmr.READYSET_ACTIONS``): dictionaries for
    unknown residues only, no hydrogens, no geometry optimisation. NASolve's
    curated dictionaries for curated residues in the model (C38, 5IU, DU, ...)
    and any ``--cif`` files are given to ReadySet, combined into
    ``input_dictionaries.cif``. ReadySet's model must keep every atom with its
    coordinates, occupancy and B (it becomes ``readyset_model.pdb`` beside
    ``directory``), and must add no hydrogen. The restraints for phenix.refine
    follow PostMR's precedence (``ligand_profiles.effective_restraints``):
    ReadySet's generated dictionary, with the curated ones where NASolve
    overrides it; a ``--cif`` file whose residues ReadySet regenerated is not
    passed twice.
    """
    from .curated_ligands import CURATED_LIGANDS, ligand_dictionary, validate_ligand_dictionary
    from .ligand_profiles import (
        combine_dictionary_inputs,
        effective_restraints,
        normalize_ccp4_torsion_alternates,
    )
    from .phosphate import PhosphateError
    from .postmr import READYSET_ACTIONS

    directory.mkdir(parents=True)
    residues = {
        line[17:20].strip()
        for line in start_model.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    }
    curated_codes = tuple(sorted(residues & set(CURATED_LIGANDS)))
    curated_copies: list[Path] = []
    for code in curated_codes:
        try:
            source = ligand_dictionary(code).resolve()
            validate_ligand_dictionary(code, source)
        except (KeyError, FileNotFoundError, ValueError) as exc:
            raise AutoFDPError(f"Curated dictionary for {code}: {exc}") from exc
        destination = directory / source.name
        if code == "1AP":
            normalize_ccp4_torsion_alternates(source, destination)
        else:
            shutil.copyfile(source, destination)
        curated_copies.append(destination)
    user_copies: list[Path] = []
    for path in extra_dictionaries:
        destination = directory / f"user_{path.name}"
        if destination.exists():
            raise AutoFDPError(f"Two --cif files are named {path.name}")
        shutil.copyfile(path, destination)
        user_copies.append(destination)
    combined: Path | None = None
    if curated_copies or user_copies:
        combined = directory / "input_dictionaries.cif"
        try:
            combine_dictionary_inputs([*curated_copies, *user_copies], combined)
        except PhosphateError as exc:
            raise AutoFDPError(f"The restraint dictionaries cannot be combined: {exc}") from exc
    command = [str(ready_set_executable), str(start_model), *READYSET_ACTIONS]
    if combined is not None:
        command.append(f"ready_set.input.cif_file_name={combined}")
    log = directory / "ready_set.log"
    try:
        completed = subprocess.run(
            command, cwd=directory,
            env=dict(environment) if environment is not None else None,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
        )
    except OSError as exc:
        raise AutoFDPError(f"Could not launch ReadySet: {exc}") from exc
    log.write_text(completed.stdout or "", encoding="utf-8")
    if completed.returncode:
        raise AutoFDPError(f"ReadySet failed with status {completed.returncode}; see {log}")
    updated = directory / f"{start_model.stem}.updated.pdb"
    generated = directory / f"{start_model.stem}.ligands.cif"
    if updated.is_file():
        source, mode = updated, "updated model"
    elif re.search(r"\bNo unknown residues\b", completed.stdout or "", re.I):
        source, mode = start_model, "no unknown residues"
    else:
        raise AutoFDPError(f"ReadySet wrote no updated model; see {log}")
    hydrogens = _hydrogen_count(source)
    if hydrogens:
        raise AutoFDPError(f"ReadySet added {hydrogens} hydrogen atom(s) despite hydrogens=False")
    changes = held_parameter_changes(start_model, source)
    if (
        changes["atoms_only_in_start"] or changes["atoms_only_in_refined"]
        or changes["max_coordinate_shift"] > 0.001
        or changes["max_occupancy_change"] > 0.001
        or changes["max_b_change"] > 0.01
    ):
        raise AutoFDPError(
            "ReadySet changed the model (atoms added, removed or renamed, or coordinates, "
            f"occupancy or B moved): {dict(changes)}; see {log}"
        )
    model = directory.parent / "readyset_model.pdb"
    shutil.copyfile(source, model)
    generated_path = generated if generated.is_file() else None
    dictionaries, _ = effective_restraints(curated_copies, generated_path, directory)
    regenerated = _component_codes(generated_path) if generated_path else set()
    for path in user_copies:
        codes = _component_codes(path)
        if codes and codes <= regenerated:
            continue   # ReadySet's dictionary already carries these residues
        if codes & regenerated:
            raise AutoFDPError(
                f"{path.name[5:]} defines {', '.join(sorted(codes))}, but ReadySet regenerated "
                f"only {', '.join(sorted(codes & regenerated))}; give those residues in "
                "separate --cif files"
            )
        dictionaries.append(path)
    return ReadySetResult(
        model=model,
        dictionaries=tuple(dictionaries),
        generated=generated_path,
        input_dictionaries=combined,
        curated_codes=curated_codes,
        log=log,
        command=tuple(command),
        mode=mode,
        changes=changes,
    )


@dataclass(frozen=True)
class FdpSiteResult:
    """One anomalous group after an f''-only round."""

    selection: str
    site: str
    atom_name: str
    element: str
    calibrated: bool              # False: refined too, but excluded by --site
    occupancy: float | None
    b_factor: float | None
    f_prime: float                # fixed at theory
    f_double_prime_start: float
    f_double_prime_refined: float | None
    f_double_prime_theory: float
    theory_source: str
    deviation: float | None       # refined - theory
    within_tolerance: bool | None
    alternate: str | None = None  # altloc of the site atom, if any


@dataclass(frozen=True)
class FdpRoundResult:
    """AutoFDP step 2: one f''-only phenix.refine round."""

    directory: Path
    cycles: int
    tolerance: float
    returncode: int
    r_work: float | None
    r_free: float | None
    sites: tuple[FdpSiteResult, ...]
    warnings: tuple[str, ...]
    status: str
    message: str
    parameters: Path
    log_path: Path
    model_path: Path | None
    map_coefficients: Path | None
    held_parameter_changes: Mapping[str, object] | None


def _newest(directory: Path, pattern: str, keep: Callable[[Path], bool] = lambda path: True) -> Path | None:
    found = sorted(
        (path.resolve() for path in directory.glob(pattern) if path.is_file() and keep(path)),
        key=lambda path: (path.stat().st_mtime_ns, path.name),
    )
    return found[-1] if found else None


def _run_phenix_refine(
    command: list[str],
    directory: Path,
    environment: Mapping[str, str] | None,
) -> int:
    """``--dry_run`` first, then the refinement; logs in ``directory``."""
    log_path = directory / "phenix.refine.log"
    preflight_log = directory / "phenix.refine.preflight.log"
    env = dict(environment) if environment is not None else None
    try:
        with preflight_log.open("w", encoding="utf-8") as handle:
            preflight = subprocess.run(
                [command[0], "--dry_run", *command[1:]], cwd=directory, env=env,
                text=True, stdout=handle, stderr=subprocess.STDOUT, check=False,
            )
        if preflight.returncode:
            log_path.write_text(
                f"phenix.refine refused the parameters (status {preflight.returncode}); "
                f"see {preflight_log.name}\n",
                encoding="utf-8",
            )
            return preflight.returncode
        with log_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(
                command, cwd=directory, env=env,
                text=True, stdout=handle, stderr=subprocess.STDOUT, check=False,
            )
        return completed.returncode
    except OSError as exc:
        log_path.write_text(f"NASolve could not launch phenix.refine: {exc}\n", encoding="utf-8")
        return 127


def fdp_round(
    start_model: Path,
    reflections: ReflectionInputs,
    refine_executable: Path,
    directory: Path,
    *,
    phenix_version: str,
    theory: Mapping[str, Mapping[str, object]],
    calibrated: Sequence[str] = (),
    cycles: int = 3,
    tolerance: float = 1.0,
    restraints: Sequence[Path] = (),
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> FdpRoundResult:
    """AutoFDP step 2: one f''-only phenix.refine round in ``directory``.

    Every strategy is off except ``group_anomalous``: no coordinate,
    real-space, B or occupancy refinement for any atom, and no weight
    optimisation. Every metal site of ``start_model`` gets an anomalous group
    with f' fixed at ``theory`` for its element and f'' refined from 0, for
    ``cycles`` macrocycles. The parameters and command come from AutoRefine's
    writers, so reflection labels, Free-R flags and the Phenix 2.x
    data-manager handling are the same as in AutoRefine. ``restraints`` are
    extra dictionaries (CIF) for residues phenix.refine does not know.

    Afterwards each refined f'' is compared with its theoretical value, and
    coordinates, occupancy and B are checked to be unchanged. Status
    ``AUTOFDP_FDP_ROUND_DONE`` when phenix.refine reported a refined f'' for
    every group, otherwise ``AUTOFDP_FAILED``.
    """
    from .autorefine import (
        AutoRefineError,
        ReflectionPlan,
        _final_anomalous_values,
        build_refine_command,
        parse_refinement_statistics,
        reflection_selector_policy,
        write_recipe_parameters,
    )

    if isinstance(cycles, bool) or not isinstance(cycles, int) or cycles < 1:
        raise AutoFDPError("--cycles must be a whole number of at least 1")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise AutoFDPError("--tolerance must be a positive number of electrons")
    sites, _ = _metal_sites(start_model, ())
    groups: dict[str, dict[str, float]] = {}
    described: dict[str, tuple[str, str, str, str | None]] = {}
    for candidate, _site in sites:
        element = str(candidate["element"]).upper()
        if element not in theory:
            raise AutoFDPError(f"No theoretical scattering for {element}")
        alternate = candidate.get("alternate") or None
        selection = _site_selection(str(candidate["site"]), str(candidate["atom_name"]), alternate)
        groups[selection] = {"f_prime": float(theory[element]["f_prime"]), "f_double_prime": 0.0}
        described[selection] = (str(candidate["site"]), str(candidate["atom_name"]), element, alternate)
    unknown = sorted(set(calibrated) - set(groups))
    if unknown:
        raise AutoFDPError("Not a metal site of the start model: " + "; ".join(unknown))

    plan = ReflectionPlan(
        observation_labels=reflections.observation_labels,
        free_r_label=reflections.free_r_label,
        free_r_test_value=reflections.free_r_test_value,
        anomalous=True,
        phase_labels=(),
        phase_file=None,
        anomalous_fallback=False,
        fallback_is_error=False,
        resolution_limit=reflections.resolution_limit,
        wavelength=reflections.wavelength,
        wavelength_source=dict(reflections.wavelength_source),
    )
    directory.mkdir(parents=True, exist_ok=True)
    parameters = directory / "autofdp.params"
    try:
        strategies = write_recipe_parameters(
            parameters,
            selector_policy=reflection_selector_policy(phenix_version),
            observations=reflections.mtz,
            reflection_plan=plan,
            macro_cycles=cycles,
            processor_count=max(1, processor_count or os.cpu_count() or 1),
            anomalous_atom_selections=tuple(groups),
            real_space_sites=False,
            adp_mode="none",
            refine_occupancies=False,
            anomalous_mode="fdp-only",
            refine_coordinates=False,
            anomalous_groups=groups,
        )
    except AutoRefineError as exc:
        raise AutoFDPError(str(exc)) from exc
    if tuple(strategies) != ("group_anomalous",):
        raise AutoFDPError("Internal error: the f'' round must refine group_anomalous only")
    command = build_refine_command(
        refine_executable, start_model, reflections.mtz,
        [path.expanduser().resolve() for path in restraints], parameters, plan,
    )
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n", encoding="utf-8")
    returncode = _run_phenix_refine(command, directory, environment)

    log_path = directory / "phenix.refine.log"
    log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
    model_path = _newest(directory, "refined_*.pdb")
    map_coefficients = _newest(directory, "refined_*map_coeffs.mtz")
    if map_coefficients is None and model_path is not None and model_path.with_suffix(".mtz").is_file():
        map_coefficients = model_path.with_suffix(".mtz")
    statistics = parse_refinement_statistics(
        log_text, model_path.read_text(encoding="utf-8", errors="replace") if model_path else ""
    )
    refined_values = _final_anomalous_values(log_text)
    warnings: list[str] = []
    changes = held_parameter_changes(start_model, model_path) if model_path is not None else None
    if changes is not None and (
        changes["max_coordinate_shift"] > 0.01
        or changes["max_occupancy_change"] > 0.01
        or changes["max_b_change"] > 0.1
    ):
        warnings.append(
            "Held parameters changed during the f'' round: "
            f"coordinates up to {changes['max_coordinate_shift']} A, occupancy up to "
            f"{changes['max_occupancy_change']}, B up to {changes['max_b_change']}"
        )
    records = _coordinate_records(start_model)
    results: list[FdpSiteResult] = []
    for selection, (site, atom_name, element, alternate) in described.items():
        record = records.get(_site_key(site, atom_name, alternate))
        refined = refined_values.get(selection)
        theory_fdp = float(theory[element]["f_double_prime"])
        deviation = None if refined is None else round(refined - theory_fdp, 3)
        results.append(FdpSiteResult(
            selection=selection,
            site=site,
            atom_name=atom_name,
            element=element,
            calibrated=not calibrated or selection in calibrated,
            occupancy=record[3] if record else None,
            b_factor=record[4] if record else None,
            f_prime=float(theory[element]["f_prime"]),
            f_double_prime_start=0.0,
            f_double_prime_refined=refined,
            f_double_prime_theory=theory_fdp,
            theory_source=str(theory[element]["source"]),
            deviation=deviation,
            within_tolerance=None if deviation is None else abs(deviation) <= tolerance,
            alternate=alternate,
        ))

    if returncode != 0 or model_path is None:
        status = STATUS_FAILED
        message = f"phenix.refine failed (exit status {returncode}); see {log_path}"
    elif any(site.f_double_prime_refined is None for site in results):
        status = STATUS_FAILED
        message = "phenix.refine reported no refined f'' for " + ", ".join(
            site.selection for site in results if site.f_double_prime_refined is None
        )
    else:
        status = STATUS_ROUND_DONE
        outside = [site for site in results if site.within_tolerance is False]
        message = (
            f"refined f'' within {tolerance:g} e- of theory for every site"
            if not outside else
            f"{len(outside)} of {len(results)} site(s) more than {tolerance:g} e- from theory"
        )
    return FdpRoundResult(
        directory=directory,
        cycles=cycles,
        tolerance=tolerance,
        returncode=returncode,
        r_work=statistics.get("r_work"),
        r_free=statistics.get("r_free"),
        sites=tuple(results),
        warnings=tuple(warnings),
        status=status,
        message=message,
        parameters=parameters,
        log_path=log_path,
        model_path=model_path,
        map_coefficients=map_coefficients,
        held_parameter_changes=changes,
    )


@dataclass(frozen=True)
class OccupancyStep:
    """The next occupancy proposed for one site (``propose_occupancy``)."""

    occupancy: float | None   # None: the site is already at the bound it must pass
    method: str               # "proportional", "secant" or "none"
    bound: str | None = None  # "floor" or "ceiling" when the step stopped there
    note: str = ""


def propose_occupancy(
    history: Sequence[tuple[float, float | None]],
    theory: float,
    *,
    floor: float = 0.05,
    ceiling: float = 1.0,
    max_ratio: float = OCCUPANCY_MAX_RATIO,
) -> OccupancyStep:
    """The next occupancy of one site, from its (occupancy, refined f'') rounds.

    The anomalous signal fixes occupancy x f'', so the refined f'' rises when
    the modelled occupancy q falls. f'' below ``theory`` means the site holds
    less metal than modelled: lower q. Above: raise q.

    - First step for the site (no earlier round with another occupancy):
      proportional, q' = q x f''ref / f''th, the occupancy at which q x f''
      would stay the same.
    - Later steps: a secant on 1/q through the site's last two rounds with
      different occupancy, fitting f'' = a + S/q and solving for f'' = theory
      (exact when f'' is inversely proportional to q up to an offset). Used
      only when f'' rose as q fell (S > 0) and the fit reaches theory at a
      positive occupancy; otherwise proportional again.
    - One step changes q at most ``max_ratio``-fold and at least 0.01 (the PDB
      occupancy column), always in the direction above, and q' is rounded to
      0.01 (half up; towards q when the step limit applies) within
      [``floor``, ``ceiling``].
    - A site already at the floor (or ceiling) that must go further gets
      ``occupancy=None``: "f'' below theory at the occupancy floor" (or "f''
      above theory at full occupancy").
    """
    if not history:
        raise AutoFDPError("No f'' round for this site yet")
    occupancy, refined = history[-1]
    if refined is None or not math.isfinite(refined):
        raise AutoFDPError("The last round gave no refined f'' for this site")
    if not (math.isfinite(theory) and theory > 0):
        raise AutoFDPError("The theoretical f'' must be a positive number")
    if not (math.isfinite(floor) and 0 < floor < ceiling <= 1.0):
        raise AutoFDPError("The occupancy floor must be above 0 and below the ceiling (at most 1)")
    if not (math.isfinite(max_ratio) and max_ratio > 1):
        raise AutoFDPError("The step limit must be more than 1-fold")
    if refined == theory:
        return OccupancyStep(occupancy, "none", note="f'' equals theory")
    lower = refined < theory
    current = round(occupancy * 100)
    lowest = math.ceil(floor * 100 - 1e-6)
    highest = math.floor(ceiling * 100 + 1e-6)
    if lower and current <= lowest:
        return OccupancyStep(None, "none", "floor", "f'' below theory at the occupancy floor")
    if not lower and current >= highest:
        reason = "f'' above theory at full occupancy" if highest == 100 else (
            "f'' above theory at the occupancy ceiling")
        return OccupancyStep(None, "none", "ceiling", reason)

    notes: list[str] = []
    method = "proportional"
    target = occupancy * refined / theory if refined > 0 else 0.0
    earlier = next(
        (point for point in reversed(history[:-1])
         if point[1] is not None and math.isfinite(point[1]) and abs(point[0] - occupancy) >= 0.005),
        None,
    )
    if earlier is not None:
        earlier_occupancy, earlier_refined = earlier
        slope = (refined - earlier_refined) / (1 / occupancy - 1 / earlier_occupancy)
        if slope <= 0:
            notes.append("secant not used: f'' did not rise as occupancy fell")
        else:
            inverse = 1 / occupancy + (theory - refined) / slope
            if inverse <= 0:
                notes.append("secant not used: the fit reaches theory at no occupancy")
            else:
                target, method = 1 / inverse, "secant"
    low, high = occupancy / max_ratio, occupancy * max_ratio
    limited = not low <= target <= high
    if limited:
        target = min(max(target, low), high)
        notes.append(f"step limited to {max_ratio:g}-fold")
        # Round towards the current occupancy so the step stays within the limit.
        proposed = math.ceil(target * 100 - 1e-6) if lower else math.floor(target * 100 + 1e-6)
    else:
        proposed = math.floor(target * 100 + 0.5)
    proposed = min(proposed, current - 1) if lower else max(proposed, current + 1)
    bound = None
    if proposed < lowest:
        proposed, bound = lowest, "floor"
        notes.append("stopped at the occupancy floor")
    elif proposed > highest:
        proposed, bound = highest, "ceiling"
        notes.append("stopped at full occupancy" if highest == 100 else "stopped at the ceiling")
    return OccupancyStep(proposed / 100, method, bound, "; ".join(notes))


@dataclass(frozen=True)
class WrongResponse:
    """A site whose f'' moved with its occupancy in its last adjustment."""

    occupancy_before: float
    occupancy_after: float
    f_double_prime_before: float
    f_double_prime_after: float
    lowered: bool                  # occupancy went down (and f'' fell)
    best_occupancy: float          # earlier occupancy to go back to
    best_f_double_prime: float     # its f'': highest (lowered) or lowest (raised)


def wrong_response(history: Sequence[tuple[float | None, float | None]]) -> WrongResponse | None:
    """Did f'' fall as occupancy fell (or rise as it rose) in the site's last round?

    The anomalous signal fixes occupancy x f'', so f'' should move against
    the occupancy. When the site's last two rounds have different occupancies
    and f'' moved with it instead, return the earlier occupancy to go back to:
    of every round of the site at another occupancy, the one with the highest
    f'' when occupancy was lowered (the lowest f'' when it was raised).
    """
    if len(history) < 2:
        return None
    (before, f_before), (after, f_after) = history[-2], history[-1]
    if None in (before, after, f_before, f_after) or abs(after - before) < 0.005:
        return None
    lowered = after < before
    if (lowered and f_after >= f_before) or (not lowered and f_after <= f_before):
        return None
    earlier = [
        (occupancy, refined) for occupancy, refined in history
        if occupancy is not None and refined is not None and abs(occupancy - after) >= 0.005
    ]
    choose = max if lowered else min
    best_occupancy, best_refined = choose(earlier, key=lambda point: point[1])
    return WrongResponse(before, after, f_before, f_after, lowered, best_occupancy, best_refined)


def _set_occupancy(line: str, value: float) -> str:
    """Replace the occupancy columns (55-60) and keep the rest of the line."""
    if not 0 <= value <= 1:
        raise AutoFDPError(f"Occupancy {value} is outside 0-1")
    body = line.rstrip("\r\n")
    ending = line[len(body):]
    body = body.ljust(60)
    return body[:54] + f"{value:6.2f}" + body[60:] + ending


def write_site_occupancy(
    source: Path,
    destination: Path,
    site: str,
    atom_name: str,
    alternate: str | None,
    occupancy: float,
) -> None:
    """Copy ``source`` with one atom's occupancy changed; nothing else."""
    key = _site_key(site, atom_name, alternate)
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    matches = []
    for index, line in enumerate(lines):
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM  ", "HETATM")) and _record_key(line) == key:
            matches.append(index)
    if len(matches) != 1:
        raise AutoFDPError(
            f"{site} {atom_name} matches {len(matches)} atoms of {source.name}; expected 1"
        )
    lines[matches[0]] = _set_occupancy(lines[matches[0]], occupancy)
    destination.write_text("".join(lines), encoding="utf-8", newline="")


@dataclass(frozen=True)
class OccupancyAdjustment:
    """One step 3 change: the occupancy one site had in round ``round``.

    ``method`` is "proportional" or "secant" for an adjustment, "restore" when
    the site went back to an earlier occupancy (``wrong_response``).
    """

    round: int
    selection: str
    site: str
    atom_name: str
    occupancy_before: float
    occupancy_after: float
    method: str
    bound: str | None
    note: str


@dataclass(frozen=True)
class SiteOutcome:
    """Where one calibrated site ended up after step 3."""

    selection: str
    site: str
    atom_name: str
    element: str
    occupancy_start: float | None
    occupancy_final: float | None
    adjustments: int
    f_double_prime_refined: float | None   # in the last round
    f_double_prime_theory: float
    deviation: float | None
    converged: bool
    reason: str
    history: tuple[tuple[int, float | None, float | None], ...]   # (round, occupancy, f'')
    held: int = 0          # times it went back to an earlier occupancy (restores)


@dataclass(frozen=True)
class AutoFDPPlan:
    """Everything checked before anything is written (the preflight)."""

    model: Path
    reflections: ReflectionInputs
    b_start: BStart
    theory: Mapping[str, Mapping[str, object]]
    selections: tuple[str, ...]          # every metal site: one anomalous group each
    calibrated: tuple[str, ...]          # the sites whose B step 1 sets
    restraints: tuple[Path, ...]
    output_root: Path                    # run_NNN is created inside it


@dataclass(frozen=True)
class AutoFDPResult:
    plan: AutoFDPPlan
    run_directory: Path
    b_start: BStart
    round: FdpRoundResult                              # the last round run
    readyset: ReadySetResult
    report_path: Path
    log_path: Path
    status: str
    message: str
    rounds: tuple[FdpRoundResult, ...] = ()            # every round, round_001 first
    adjustments: tuple[OccupancyAdjustment, ...] = ()  # step 3, in order
    outcomes: tuple[SiteOutcome, ...] = ()             # each calibrated site at the end
    final_model: Path | None = None                    # the last round's model


def plan_autofdp(
    model: Path,
    mtz: Path,
    refine_executable: Path,
    mtz_dump_executable: Path,
    *,
    sites: Sequence[str] = (),
    b_start: str | float = "surrounding",
    theory_overrides: Mapping[str, float] | None = None,
    wavelength: float | None = None,
    restraints: Sequence[Path] = (),
    output_root: Path | None = None,
    environment: Mapping[str, str] | None = None,
) -> AutoFDPPlan:
    """Preflight: check the inputs and return the plan, writing nothing.

    Reads the model's metal sites and their planned B (step 1), the MTZ
    labels, Free-R flags and wavelength (``phenix.mtz.dump``), and the
    theoretical f'/f'' of every metal element (``phenix.python``). Runs no
    refinement. ``output_root`` defaults to ``AutoFDP`` beside the MTZ file.
    """
    model = checked_model(model)
    plan_b = prepare_start_model(model, sites=sites, b_start=b_start)
    found, _ = _metal_sites(model, ())
    selections = tuple(
        _site_selection(str(candidate["site"]), str(candidate["atom_name"]),
                        candidate.get("alternate") or None)
        for candidate, _ in found
    )
    calibrated = tuple(
        _site_selection(site.site, site.atom_name, site.alternate) for site in plan_b.sites
    )
    files = []
    for path in restraints:
        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            raise AutoFDPError(f"No restraint file: {resolved}")
        files.append(resolved)
    reflections = read_reflections(mtz, mtz_dump_executable, environment, wavelength=wavelength)
    elements = {str(candidate["element"]).upper() for candidate, _ in found}
    theory = theoretical_scattering(
        sorted(elements), reflections.wavelength, refine_executable, environment, theory_overrides
    )
    root = (output_root or reflections.mtz.parent / "AutoFDP").expanduser().resolve()
    return AutoFDPPlan(
        model=model,
        reflections=reflections,
        b_start=plan_b,
        theory=theory,
        selections=selections,
        calibrated=calibrated,
        restraints=tuple(files),
        output_root=root,
    )


def _next_run_directory(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    numbers = [
        int(match.group(1))
        for path in root.glob("run_*")
        if path.is_dir() and (match := re.fullmatch(r"run_(\d+)", path.name))
    ]
    destination = root / f"run_{max(numbers, default=0) + 1:03d}"
    destination.mkdir()
    return destination


def _file_record(path: Path, run_directory: Path) -> dict[str, object]:
    """A portable reference (run or dataset anchor) with checksum and size.

    The absolute path at creation is kept as provenance only.
    """
    return {
        **artifact_reference(path, run_directory),
        "absolute_path_at_creation": str(path.resolve()),
        "sha256": file_sha256(path),
        "size": path.stat().st_size,
    }


def _site_lines(round_result: FdpRoundResult) -> list[str]:
    lines = []
    for site in round_result.sites:
        refined = "none" if site.f_double_prime_refined is None else f"{site.f_double_prime_refined:.2f}"
        deviation = "n/a" if site.deviation is None else f"{site.deviation:+.2f}"
        occupancy = "-" if site.occupancy is None else f"{site.occupancy:.2f}"
        lines.append(
            f"{site.site} {site.atom_name} ({site.element}): occupancy {occupancy}; "
            f"f' {site.f_prime:.2f} fixed; f'' {refined} vs theory "
            f"{site.f_double_prime_theory:.2f} ({site.theory_source}); deviation {deviation}"
        )
    return lines


def _log_lines(
    plan: AutoFDPPlan,
    b_start: BStart,
    rounds: Sequence[FdpRoundResult],
    readyset: "ReadySetResult | None" = None,
    adjustments: Sequence[OccupancyAdjustment] = (),
    outcomes: Sequence[SiteOutcome] = (),
) -> list[str]:
    reflections = plan.reflections
    lines = [
        f"AutoFDP, {_now()}",
        f"Model: {plan.model}",
        f"MTZ: {reflections.mtz}",
        f"Observations: {','.join(reflections.observation_labels)}; Free-R "
        f"{reflections.free_r_label} (test value {reflections.free_r_test_value})",
        f"Wavelength: {reflections.wavelength:g} A ({reflections.wavelength_source.get('method')})",
        "",
        "Step 1 (B start)",
    ]
    for site in b_start.sites:
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
    lines.extend(f"Warning: {warning}" for warning in b_start.warnings)
    if readyset is not None:
        lines.extend([
            "",
            f"ReadySet (no hydrogens): {readyset.mode}; model unchanged "
            f"({readyset.changes.get('atoms_compared')} atoms compared)",
            "Curated dictionaries: " + (", ".join(readyset.curated_codes) or "none"),
            "Refinement dictionaries: "
            + (", ".join(path.name for path in readyset.dictionaries) or "none"),
        ])
    if rounds:
        first = rounds[0]
        lines.extend(["", f"Step 2 (f''-only round, round_001): {first.status}: {first.message}",
                      f"R-work {first.r_work}, R-free {first.r_free}"])
        lines.extend(_site_lines(first))
        lines.extend(f"Warning: {warning}" for warning in first.warnings)
    if len(rounds) > 1 or outcomes:
        lines.extend(["", "Step 3 (occupancy adjustment; same f''-only round each time)"])
        changed = {adjustment.round: adjustment for adjustment in adjustments}
        for number, result in enumerate(rounds[1:], start=2):
            adjustment = changed.get(number)
            what = "" if adjustment is None else (
                f"{adjustment.site} {adjustment.atom_name} occupancy "
                f"{adjustment.occupancy_before:.2f} -> {adjustment.occupancy_after:.2f} "
                f"({adjustment.method}{'; ' + adjustment.note if adjustment.note else ''}); "
            )
            lines.append(f"round_{number:03d}: {what}{result.status}; "
                         f"R-work {result.r_work}, R-free {result.r_free}")
            lines.extend("  " + line for line in _site_lines(result))
            lines.extend(f"  Warning: {warning}" for warning in result.warnings)
            if result.status != STATUS_ROUND_DONE:
                lines.append(f"  {result.message}")
        for outcome in outcomes:
            start = "-" if outcome.occupancy_start is None else f"{outcome.occupancy_start:.2f}"
            final = "-" if outcome.occupancy_final is None else f"{outcome.occupancy_final:.2f}"
            refined = ("none" if outcome.f_double_prime_refined is None
                       else f"{outcome.f_double_prime_refined:.2f}")
            lines.append(
                f"{outcome.site} {outcome.atom_name} ({outcome.element}): occupancy {start} -> "
                f"{final}; f'' {refined} vs theory {outcome.f_double_prime_theory:.2f}; "
                f"{outcome.reason}"
            )
    return lines


def _round_entry(
    number: int,
    round_result: FdpRoundResult,
    input_model: Path,
    run_directory: Path,
    theory: Mapping[str, Mapping[str, object]],
    cycles: int,
    adjustment: OccupancyAdjustment | None = None,
) -> dict[str, object]:
    return {
        "round": number,
        "directory": f"round_{number:03d}",
        "strategy": ["group_anomalous"],
        "input_model": _file_record(input_model, run_directory),
        "adjusted_site": None if adjustment is None else adjustment.selection,
        "occupancy_before": None if adjustment is None else adjustment.occupancy_before,
        "occupancy_after": None if adjustment is None else adjustment.occupancy_after,
        "method": None if adjustment is None else adjustment.method,
        "note": None if adjustment is None else adjustment.note,
        "starting_groups": {
            site.selection: {"f_prime": float(theory[site.element]["f_prime"]),
                             "f_double_prime": 0.0}
            for site in round_result.sites
        },
        "cycles": cycles,
        "returncode": round_result.returncode,
        "r_work": round_result.r_work,
        "r_free": round_result.r_free,
        "held_parameter_changes": round_result.held_parameter_changes,
        "refined_model": (
            _file_record(round_result.model_path, run_directory) if round_result.model_path else None
        ),
        "map_coefficients": (
            _file_record(round_result.map_coefficients, run_directory)
            if round_result.map_coefficients else None
        ),
        "sites": [asdict(site) for site in round_result.sites],
        "warnings": list(round_result.warnings),
        "status": round_result.status,
        "message": round_result.message,
    }


def _site_outcomes(
    rounds: Sequence[FdpRoundResult],
    adjustments: Sequence[OccupancyAdjustment],
    stopped: Mapping[str, str],
    tolerance: float,
) -> list[SiteOutcome]:
    """Each calibrated site in the last round, with its history and why it stopped."""
    last = rounds[-1]
    outcomes = []
    for site in last.sites:
        if not site.calibrated:
            continue
        history = []
        for number, result in enumerate(rounds, start=1):
            match = next((item for item in result.sites if item.selection == site.selection), None)
            if match is not None:
                history.append((number, match.occupancy, match.f_double_prime_refined))
        count = sum(1 for adjustment in adjustments
                    if adjustment.selection == site.selection and adjustment.method != "restore")
        held = sum(1 for adjustment in adjustments
                   if adjustment.selection == site.selection and adjustment.method == "restore")
        converged = bool(site.within_tolerance)
        if converged:
            reason = (f"within +/-{tolerance:g} e- of theory"
                      + (f" after {count} adjustment(s)" if count else " without adjustment")
                      + (f", held back {held} time(s)" if held else ""))
        elif site.selection in stopped and not stopped[site.selection].startswith("within"):
            reason = stopped[site.selection]
        elif site.selection in stopped:
            reason = (f"was within +/-{tolerance:g} e-, then moved to "
                      f"{site.deviation:+.2f} while other sites were adjusted")
        else:
            reason = "not adjusted"
        outcomes.append(SiteOutcome(
            selection=site.selection,
            site=site.site,
            atom_name=site.atom_name,
            element=site.element,
            occupancy_start=history[0][1] if history else None,
            occupancy_final=site.occupancy,
            adjustments=count,
            f_double_prime_refined=site.f_double_prime_refined,
            f_double_prime_theory=site.f_double_prime_theory,
            deviation=site.deviation,
            converged=converged,
            reason=reason,
            history=tuple(history),
            held=held,
        ))
    return outcomes


def execute_autofdp(
    plan: AutoFDPPlan,
    refine_executable: Path,
    ready_set_executable: Path,
    *,
    phenix_version: str,
    b_start: str | float = "surrounding",
    sites: Sequence[str] = (),
    cycles: int = 3,
    tolerance: float = 1.0,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    occupancy_floor: float = 0.05,
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> AutoFDPResult:
    """Run steps 1 to 3 in a new ``run_NNN`` folder of ``plan.output_root``.

    Step 1 writes ``source_model.pdb`` and ``start_model.pdb``. Step 2 runs
    ReadySet on the start model (``ReadySet/``, no hydrogens) for the restraint
    dictionaries, then the f''-only round on ReadySet's model in
    ``round_001``.

    Step 3 adjusts occupancy, one calibrated site at a time: the site whose
    refined f'' is furthest from theory (outside ``tolerance``) first. Its
    occupancy is changed (``propose_occupancy``; lower when f'' is below
    theory, higher when above, within [``occupancy_floor``, 1]) in a copy of
    the previous round's model (``round_NNN/input_model.pdb``; nothing else
    changes), and the same f''-only round is run again on it. This repeats
    until the site is within ``tolerance``, has had ``max_iterations``
    adjustments, or is at the occupancy bound it would have to pass; then the
    next site, ranked again by the latest round. Every round refines every
    site's f'' from 0 with the same parameters. ``max_iterations=0`` runs
    round_001 only.

    When an adjustment moves f'' with the occupancy (``wrong_response``:
    f'' fell as occupancy fell, or rose as it rose), the site goes back to its
    earlier occupancy with the highest (lowest) f'' in a "restore" round and
    is held there while the other sites are adjusted. When no other site can
    be adjusted, it is taken up again (same per-site limit). A second such
    response stops it at that occupancy, again with a restore round.

    Status ``AUTOFDP_CONVERGED`` when every calibrated site is within
    ``tolerance`` in the last round (written as ``final_model.pdb``),
    ``AUTOFDP_REVIEW`` when one is not, ``AUTOFDP_FAILED`` when a round
    fails. ``report.json`` (workflow ``autofdp``, schema 2) and
    ``autofdp.log`` record the inputs (paths and checksums), what was read
    from the MTZ, the theory, the B changes, every round and the occupancy
    history of every site. Nothing outside the new folder is written. The
    identity check, B confirmation and final round are not implemented.
    """
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 0:
        raise AutoFDPError("--max-iterations must be a whole number of at least 0")
    if not (math.isfinite(occupancy_floor) and 0 < occupancy_floor < 1):
        raise AutoFDPError("--occupancy-floor must be above 0 and below 1")
    run_directory = _next_run_directory(plan.output_root)
    written = prepare_start_model(
        plan.model, sites=sites, b_start=b_start, output_directory=run_directory
    )
    assert written.start_model is not None and written.source_copy is not None
    reflections = plan.reflections
    report: dict[str, object] = {
        "workflow": WORKFLOW,
        "schema_version": SCHEMA_VERSION,
        "status": STATUS_PREPARED,
        "created_utc": _now(),
        "phenix_version": phenix_version,
        "inputs": {
            "model": _file_record(plan.model, run_directory),
            "mtz": _file_record(reflections.mtz, run_directory),
            "restraints": [_file_record(path, run_directory) for path in plan.restraints],
        },
        "reflections": {
            "observation_labels": list(reflections.observation_labels),
            "free_r_label": reflections.free_r_label,
            "free_r_test_value": reflections.free_r_test_value,
            "wavelength": reflections.wavelength,
            "wavelength_source": dict(reflections.wavelength_source),
            "resolution_limit": reflections.resolution_limit,
        },
        "settings": {
            "sites": list(sites), "b_start": b_start, "cycles": cycles, "tolerance": tolerance,
            "max_iterations": max_iterations, "occupancy_floor": occupancy_floor,
        },
        "theory": {element: dict(values) for element, values in plan.theory.items()},
        "b_start": {
            "source_model": _file_record(written.source_copy, run_directory),
            "start_model": _file_record(written.start_model, run_directory),
            "sites": [_site_record(site) for site in written.sites],
            "warnings": list(written.warnings),
        },
        "rounds": [],
        "steps": {"b_start": "done", "fdp_rounds": "running", "occupancy_adjustment": "pending"},
    }
    report_path = run_directory / "report.json"
    _write_json(report_path, report)
    log_path = run_directory / "autofdp.log"
    readyset: ReadySetResult | None = None
    try:
        readyset = run_readyset(
            written.start_model, run_directory / "ReadySet", ready_set_executable,
            extra_dictionaries=plan.restraints, environment=environment,
        )
        report["readyset"] = {
            "command": list(readyset.command),
            "mode": readyset.mode,
            "log": _file_record(readyset.log, run_directory),
            "model": _file_record(readyset.model, run_directory),
            "generated_dictionary": (
                _file_record(readyset.generated, run_directory) if readyset.generated else None
            ),
            "input_dictionaries": (
                _file_record(readyset.input_dictionaries, run_directory)
                if readyset.input_dictionaries else None
            ),
            "curated_residues": list(readyset.curated_codes),
            "refinement_dictionaries": [
                _file_record(path, run_directory) for path in readyset.dictionaries
            ],
            "model_changes": dict(readyset.changes),
        }
        report["steps"]["readyset"] = "done"
        _write_json(report_path, report)
    except AutoFDPError as exc:
        report.update(status=STATUS_FAILED, message=str(exc))
        report["steps"]["readyset"] = "failed"
        report["steps"]["fdp_rounds"] = "failed before phenix.refine ran"
        report["steps"]["occupancy_adjustment"] = "not run"
        _write_json(report_path, report)
        log_path.write_text("\n".join(_log_lines(plan, written, (), readyset)
                                      + ["", f"Failed: {exc}"]) + "\n", encoding="utf-8")
        raise

    rounds: list[FdpRoundResult] = []
    inputs: list[Path] = []
    adjustments: list[OccupancyAdjustment] = []

    def refine(number: int, model: Path, adjustment: OccupancyAdjustment | None) -> FdpRoundResult:
        result = fdp_round(
            model,
            reflections,
            refine_executable,
            run_directory / f"round_{number:03d}",
            phenix_version=phenix_version,
            theory=plan.theory,
            calibrated=plan.calibrated,
            cycles=cycles,
            tolerance=tolerance,
            restraints=readyset.dictionaries,
            environment=environment,
            processor_count=processor_count,
        )
        rounds.append(result)
        inputs.append(model)
        report["rounds"].append(_round_entry(
            number, result, model, run_directory, plan.theory, cycles, adjustment
        ))
        report["steps"]["fdp_rounds"] = f"{len(rounds)} run"
        report.update(status=result.status, message=f"round_{number:03d}: {result.message}")
        _write_json(report_path, report)
        return result

    try:
        first = refine(1, readyset.model, None)
    except AutoFDPError as exc:
        report.update(status=STATUS_FAILED, message=str(exc))
        report["steps"]["fdp_rounds"] = "failed before phenix.refine ran"
        report["steps"]["occupancy_adjustment"] = "not run"
        _write_json(report_path, report)
        log_path.write_text("\n".join(_log_lines(plan, written, (), readyset)
                                      + ["", f"Failed: {exc}"]) + "\n", encoding="utf-8")
        raise

    # Step 3: one site at a time, the one furthest from theory first.
    failure = None if first.status == STATUS_ROUND_DONE else first.message
    history: dict[str, list[tuple[float | None, float | None]]] = {}
    for site in first.sites:
        history[site.selection] = [(site.occupancy, site.f_double_prime_refined)]
    stopped: dict[str, str] = {}    # selection -> why its adjustment ended
    used: dict[str, int] = {}       # adjustments per site (restores not counted)
    held: dict[str, int] = {}       # times each site went back to an earlier occupancy
    holding: set[str] = set()       # held back now, while the other sites are adjusted
    active: str | None = None

    def change(site: FdpSiteResult, occupancy: float, method: str,
               bound: str | None = None, note: str = "") -> None:
        """One round with ``site`` at ``occupancy``; nothing else changes."""
        nonlocal failure
        number = len(rounds) + 1
        directory = run_directory / f"round_{number:03d}"
        directory.mkdir()
        model = directory / "input_model.pdb"
        write_site_occupancy(
            inputs[-1], model, site.site, site.atom_name, site.alternate, occupancy
        )
        adjustment = OccupancyAdjustment(
            round=number,
            selection=site.selection,
            site=site.site,
            atom_name=site.atom_name,
            occupancy_before=site.occupancy,
            occupancy_after=occupancy,
            method=method,
            bound=bound,
            note=note,
        )
        adjustments.append(adjustment)
        try:
            result = refine(number, model, adjustment)
        except AutoFDPError as exc:
            failure = str(exc)
            return
        if result.status != STATUS_ROUND_DONE:
            failure = f"round_{number:03d}: {result.message}"
            return
        for item in result.sites:
            history[item.selection].append((item.occupancy, item.f_double_prime_refined))

    while failure is None:
        latest = {site.selection: site for site in rounds[-1].sites}
        resumed = False
        if active is not None:
            site = latest[active]
            response = wrong_response(history[active])
            if site.within_tolerance:
                stopped[active] = f"within +/-{tolerance:g} e- of theory"
                active = None
            elif response is not None:
                # f'' moved with occupancy: back to the earlier occupancy with
                # the best f''. The first time, hold the site there while the
                # other sites are adjusted; the second time, stop it there.
                held[active] = held.get(active, 0) + 1
                fell = "fell" if response.lowered else "rose"
                best = "highest" if response.lowered else "lowest"
                moved = (f"{response.occupancy_before:.2f} -> {response.occupancy_after:.2f}")
                if held[active] == 1:
                    holding.add(active)
                    note = (
                        f"f'' {fell} from {response.f_double_prime_before:.2f} to "
                        f"{response.f_double_prime_after:.2f} as occupancy went {moved}; back to "
                        f"{response.best_occupancy:.2f} ({best} f'' "
                        f"{response.best_f_double_prime:.2f}), held while the other sites are adjusted"
                    )
                else:
                    stopped[active] = note = (
                        f"f'' {fell} again when occupancy was "
                        f"{'lowered' if response.lowered else 'raised'} ({moved}); kept at "
                        f"{response.best_occupancy:.2f}, where f'' was {best} "
                        f"({response.best_f_double_prime:.2f})"
                    )
                active = None
                change(site, response.best_occupancy, "restore", note=note)
                continue
            elif used[active] >= max_iterations:
                stopped[active] = (
                    f"still {site.deviation:+.2f} e- from theory after "
                    f"{max_iterations} adjustment(s) (--max-iterations)"
                )
                active = None
        if active is None:
            pending = [
                site for site in rounds[-1].sites
                if site.calibrated and site.selection not in stopped
                and site.selection not in holding and not site.within_tolerance
            ]
            if pending:
                if max_iterations == 0:
                    for site in pending:
                        stopped[site.selection] = "not adjusted (--max-iterations 0)"
                    break
                active = max(pending, key=lambda item: abs(item.deviation)).selection
                used.setdefault(active, 0)
            elif holding:
                # Nothing else to adjust: take up a held site again.
                active = max(
                    (latest[selection] for selection in holding),
                    key=lambda item: abs(item.deviation),
                ).selection
                holding.discard(active)
                resumed = True
                if latest[active].within_tolerance:
                    stopped[active] = f"within +/-{tolerance:g} e- of theory"
                    active = None
                    continue
                if used[active] >= max_iterations:
                    stopped[active] = (
                        f"still {latest[active].deviation:+.2f} e- from theory after "
                        f"{max_iterations} adjustment(s) (--max-iterations)"
                    )
                    active = None
                    continue
            else:
                break
        site = latest[active]
        if site.occupancy is None:
            failure = f"{site.site} {site.atom_name} has no readable occupancy"
            break
        step = propose_occupancy(
            history[active], site.f_double_prime_theory, floor=occupancy_floor
        )
        if step.occupancy is None:
            stopped[active] = step.note
            active = None
            continue
        note = step.note
        if resumed:
            note = "taken up again after the other sites" + (f"; {note}" if note else "")
        used[active] += 1
        change(site, step.occupancy, step.method, step.bound, note)

    finished = [result for result in rounds if result.status == STATUS_ROUND_DONE]
    outcomes = _site_outcomes(finished, adjustments, stopped, tolerance) if finished else []
    final_model: Path | None = None
    if failure is not None:
        status, message = STATUS_FAILED, failure
        report["steps"]["occupancy_adjustment"] = "failed" if adjustments else "not run"
    else:
        final_model = run_directory / "final_model.pdb"
        shutil.copyfile(inputs[-1], final_model)
        report["final_model"] = _file_record(final_model, run_directory)
        open_sites = [outcome for outcome in outcomes if not outcome.converged]
        if not open_sites:
            status = STATUS_CONVERGED
            message = (
                f"every site within +/-{tolerance:g} e- of theory after "
                f"{len(adjustments)} occupancy adjustment(s) in {len(rounds)} f'' round(s)"
            )
            report["steps"]["occupancy_adjustment"] = "done"
        else:
            status = STATUS_REVIEW
            message = "; ".join(
                f"{outcome.site} {outcome.atom_name}: {outcome.reason}" for outcome in open_sites
            )
            report["steps"]["occupancy_adjustment"] = "review"
    report["occupancy"] = {
        "method": (
            "one site at a time, furthest from theory first; proportional first step "
            "(q x f''/theory), then a secant on 1/q through the site's last two "
            "occupancies; at most 4-fold and at least 0.01 per step. When f'' moves "
            "with occupancy (falls as it falls, or rises as it rises), the site goes "
            "back to its earlier occupancy with the highest (lowest) f'' and is held "
            "there while the other sites are adjusted, then taken up again; a second "
            "such response stops it there"
        ),
        "max_iterations": max_iterations,
        "floor": occupancy_floor,
        "ceiling": 1.0,
        "max_ratio": OCCUPANCY_MAX_RATIO,
        "adjustments": [asdict(adjustment) for adjustment in adjustments],
        "sites": [
            {
                **{key: value for key, value in asdict(outcome).items() if key != "history"},
                "history": [
                    {"round": number, "occupancy": occupancy, "f_double_prime": refined,
                     "deviation": (None if refined is None
                                   else round(refined - outcome.f_double_prime_theory, 3))}
                    for number, occupancy, refined in outcome.history
                ],
            }
            for outcome in outcomes
        ],
    }
    report.update(status=status, message=message)
    _write_json(report_path, report)
    log_path.write_text(
        "\n".join(_log_lines(plan, written, rounds, readyset, adjustments, outcomes)
                  + ["", f"{status}: {message}"]) + "\n",
        encoding="utf-8",
    )
    return AutoFDPResult(
        plan=plan,
        run_directory=run_directory,
        b_start=written,
        round=rounds[-1],
        readyset=readyset,
        report_path=report_path,
        log_path=log_path,
        status=status,
        message=message,
        rounds=tuple(rounds),
        adjustments=tuple(adjustments),
        outcomes=tuple(outcomes),
        final_model=final_model,
    )


def map_value_at_site(
    map_coefficients: Path,
    labels: str,
    coordinates: tuple[float, float, float],
    phenix_python: Path,
    environment: Mapping[str, str] | None = None,
) -> float:
    """Value of a map (mFo-DFc or anomalous) at a point, in sigma units."""
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
    "ANOMALOUS_LABEL_SETS",
    "AutoFDPError",
    "AutoFDPPlan",
    "AutoFDPResult",
    "AutoFDPSettings",
    "DEFAULT_MAX_ITERATIONS",
    "BStart",
    "BStartSite",
    "FdpRoundResult",
    "FdpSiteResult",
    "Neighbour",
    "NeighbourB",
    "NOT_IMPLEMENTED_MESSAGE",
    "OCCUPANCY_MAX_RATIO",
    "OccupancyAdjustment",
    "OccupancyStep",
    "ReadySetResult",
    "ReflectionInputs",
    "SiteState",
    "SiteOutcome",
    "STATUS_CONVERGED",
    "STATUS_FAILED",
    "STATUS_PREPARED",
    "STATUS_REVIEW",
    "STATUS_ROUND_DONE",
    "WORKFLOW",
    "checked_model",
    "confirm_b_round",
    "execute_autofdp",
    "fdp_round",
    "final_round",
    "held_parameter_changes",
    "identity_check",
    "map_value_at_site",
    "parse_b_start",
    "parse_theory_overrides",
    "plan_autofdp",
    "prepare_start_model",
    "propose_occupancy",
    "read_reflections",
    "run_readyset",
    "surrounding_b",
    "theoretical_scattering",
    "write_site_occupancy",
    "wrong_response",
    "WrongResponse",
]
