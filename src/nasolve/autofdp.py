"""AutoFDP: measure anomalous-scatterer occupancy by f'' refinement.

SPECIFIED, NOT IMPLEMENTED. Every function below documents what it is
supposed to do and raises NotImplementedError. The full contract, defaults
and open questions are in docs/autofdp.md.

Method (supplements of the mmDNA papers, after Liu & Hendrickson 2013):
hold coordinates, B and occupancy; fix f' at its theoretical value for the
collection wavelength; refine f'' alone; adjust occupancy until refined f''
is within +/-1 e- of the theoretical f'' of the metal.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

WORKFLOW = "autofdp"
SCHEMA_VERSION = 1

STATUS_CONVERGED = "AUTOFDP_CONVERGED"  # exit 0
STATUS_REVIEW = "AUTOFDP_REVIEW"        # exit 2
STATUS_FAILED = "AUTOFDP_FAILED"        # exit 1

RECIPE_FDP = "AutoFDP/fdp"
RECIPE_CONFIRM_B = "AutoFDP/confirm-b"
RECIPE_FINAL = "AutoFDP/final"

NOT_IMPLEMENTED_MESSAGE = (
    "nasolve autofdp is specified but not implemented yet; nothing was run or "
    "written. See docs/autofdp.md for what it will do."
)


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
    2. ``starting_checkpoint``: occupancy capped at 1, B from
       ``settings.b_start``.
    3. Repeat up to ``settings.max_iterations``: ``fdp_round``, then
       ``decide_site_update`` for every unconverged site, and write the new
       occupancy/B as an ``autofdp-edit`` checkpoint.
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


def surrounding_b(model_path: Path, site: str, atom_name: str) -> float:
    """Mean B of the nucleotides coordinating a site.

    Coordinating nucleotides have a base N, O or S within 3.6 A of the site
    (the NARestraints metal-pair rule); if there are none, nucleotides with
    any atom within 6 A are used. No nucleotide at all is an AutoFDPError.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


def starting_checkpoint(
    run_directory: Path,
    from_checkpoint: str | None,
    sites: Sequence[SiteState],
) -> str:
    """Write the starting model as an ``autofdp-edit`` checkpoint.

    Only the site atoms change: occupancy capped at 1.0, B set as planned.
    Returns the checkpoint ID.
    """
    raise NotImplementedError(NOT_IMPLEMENTED_MESSAGE)


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
    "AutoFDPSettings",
    "SiteState",
    "NOT_IMPLEMENTED_MESSAGE",
    "STATUS_CONVERGED",
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
    "plan_autofdp",
    "starting_checkpoint",
    "surrounding_b",
    "theoretical_scattering",
]
