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
- Step 2 (``fdp_round``): one phenix.refine round on that copy with every
  strategy off except ``group_anomalous`` (f' fixed at theory, f'' refined
  from 0), then each refined f'' compared with its theoretical value.

``plan_autofdp`` is the preflight; ``execute_autofdp`` runs steps 1 and 2.
Everything after that (occupancy adjustment, identity check, B confirmation,
final round) is a documented stub that raises NotImplementedError. The full
contract, defaults and open questions are in docs/autofdp.md.

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
STATUS_ROUND_DONE = "AUTOFDP_FDP_ROUND_DONE"    # exit 0; one f''-only round
STATUS_CONVERGED = "AUTOFDP_CONVERGED"          # exit 0 (later steps)
STATUS_REVIEW = "AUTOFDP_REVIEW"                # exit 2 (later steps)
STATUS_FAILED = "AUTOFDP_FAILED"                # exit 1

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


def _coordinate_records(path: Path) -> dict[tuple[str, ...], tuple[float, ...]]:
    """(chain, resseq, icode, name, altloc) -> (x, y, z, occupancy, B)."""
    records: dict[tuple[str, ...], tuple[float, ...]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(),
               line[12:16].strip(), line[16:17].strip())
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
    round: FdpRoundResult
    report_path: Path
    log_path: Path
    status: str
    message: str


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


def _log_lines(plan: AutoFDPPlan, b_start: BStart, round_result: FdpRoundResult | None) -> list[str]:
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
    if round_result is not None:
        lines.extend(["", f"Step 2 (f''-only round): {round_result.status}: {round_result.message}",
                      f"R-work {round_result.r_work}, R-free {round_result.r_free}"])
        for site in round_result.sites:
            refined = "none" if site.f_double_prime_refined is None else f"{site.f_double_prime_refined:.2f}"
            deviation = "n/a" if site.deviation is None else f"{site.deviation:+.2f}"
            lines.append(
                f"{site.site} {site.atom_name} ({site.element}): f' {site.f_prime:.2f} fixed; "
                f"f'' {refined} vs theory {site.f_double_prime_theory:.2f} "
                f"({site.theory_source}); deviation {deviation}"
            )
        lines.extend(f"Warning: {warning}" for warning in round_result.warnings)
    return lines


def execute_autofdp(
    plan: AutoFDPPlan,
    refine_executable: Path,
    *,
    phenix_version: str,
    b_start: str | float = "surrounding",
    sites: Sequence[str] = (),
    cycles: int = 3,
    tolerance: float = 1.0,
    environment: Mapping[str, str] | None = None,
    processor_count: int | None = None,
) -> AutoFDPResult:
    """Run steps 1 and 2 in a new ``run_NNN`` folder of ``plan.output_root``.

    Step 1 writes ``source_model.pdb`` and ``start_model.pdb``; step 2 runs the
    f''-only round in ``round_001``. ``report.json`` (workflow ``autofdp``,
    schema 2) and ``autofdp.log`` record the inputs (paths and checksums),
    what was read from the MTZ, the theory, the B changes and the round.
    Nothing outside the new folder is written. The later steps (occupancy
    adjustment onward) are not implemented.
    """
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
        },
        "theory": {element: dict(values) for element, values in plan.theory.items()},
        "b_start": {
            "source_model": _file_record(written.source_copy, run_directory),
            "start_model": _file_record(written.start_model, run_directory),
            "sites": [_site_record(site) for site in written.sites],
            "warnings": list(written.warnings),
        },
        "rounds": [],
        "steps": {"b_start": "done", "fdp_rounds": "running", "occupancy_adjustment": "not implemented"},
    }
    report_path = run_directory / "report.json"
    _write_json(report_path, report)
    log_path = run_directory / "autofdp.log"
    try:
        round_result = fdp_round(
            written.start_model,
            reflections,
            refine_executable,
            run_directory / "round_001",
            phenix_version=phenix_version,
            theory=plan.theory,
            calibrated=plan.calibrated,
            cycles=cycles,
            tolerance=tolerance,
            restraints=plan.restraints,
            environment=environment,
            processor_count=processor_count,
        )
    except AutoFDPError as exc:
        report.update(status=STATUS_FAILED, message=str(exc))
        report["steps"]["fdp_rounds"] = "failed before phenix.refine ran"
        _write_json(report_path, report)
        log_path.write_text("\n".join(_log_lines(plan, written, None)
                                      + ["", f"Step 2 failed: {exc}"]) + "\n", encoding="utf-8")
        raise
    entry = {
        "round": 1,
        "directory": "round_001",
        "strategy": ["group_anomalous"],
        "starting_groups": {
            selection: {"f_prime": float(plan.theory[element]["f_prime"]), "f_double_prime": 0.0}
            for selection, element in (
                (site.selection, site.element) for site in round_result.sites
            )
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
    report["rounds"] = [entry]
    report.update(status=round_result.status, message=round_result.message)
    report["steps"]["fdp_rounds"] = "1 run"
    _write_json(report_path, report)
    log_path.write_text("\n".join(_log_lines(plan, written, round_result)) + "\n", encoding="utf-8")
    return AutoFDPResult(
        plan=plan,
        run_directory=run_directory,
        b_start=written,
        round=round_result,
        report_path=report_path,
        log_path=log_path,
        status=round_result.status,
        message=round_result.message,
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
    "ANOMALOUS_LABEL_SETS",
    "AutoFDPError",
    "AutoFDPPlan",
    "AutoFDPResult",
    "AutoFDPSettings",
    "BStart",
    "BStartSite",
    "FdpRoundResult",
    "FdpSiteResult",
    "Neighbour",
    "NeighbourB",
    "NOT_IMPLEMENTED_MESSAGE",
    "ReflectionInputs",
    "SiteState",
    "STATUS_CONVERGED",
    "STATUS_FAILED",
    "STATUS_PREPARED",
    "STATUS_REVIEW",
    "STATUS_ROUND_DONE",
    "WORKFLOW",
    "checked_model",
    "confirm_b_round",
    "decide_site_update",
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
    "read_reflections",
    "surrounding_b",
    "theoretical_scattering",
]
