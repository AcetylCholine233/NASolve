# AutoFDP: occupancy from f″ refinement

**Status: steps 1 and 2 implemented, later steps specified only (2026-11-13).**
`nasolve autofdp MODEL.pdb DATA.mtz` needs only a starting model and an MTZ
file; it reads no NASolve run report, AutoRefine round or checkpoint. Step 1
(B start) copies the model and sets each metal's B to the mean B of its
neighbouring nucleotides. Step 2 (f″ round) runs one phenix.refine round on
that copy with every strategy off except the anomalous scatterers, opens the
result and its maps in Coot, and prints each refined f″ against its
theoretical value. Occupancy adjustment and the steps after it (3-6) are
documented stubs in `src/nasolve/autofdp.py` that raise `NotImplementedError`;
options that belong to them are accepted but have no effect yet, and
`--element` is refused. The defaults marked *(to confirm)* follow my reading of
the papers below and have not yet been checked by the lab.

## Purpose

In a metal-mediated base pair, the anomalous signal of a metal site fixes the
product occupancy × f″ (reduced by the B factor), not the two separately.
`nasolve autofdp` uses that to measure occupancy and check identity:

1. Hold coordinates, B and occupancy, fix f′ at its theoretical value for the
   collection wavelength, and refine f″ alone.
2. Compare the refined f″ with the theoretical f″ of the metal at that
   wavelength (about 2.0 e⁻ for Ag and 10 e⁻ for Hg near 1.0074 Å).
3. If f″ is too high, the site holds more metal than modelled: raise the
   occupancy. If too low, lower it. Repeat until refined f″ is within ±1 e⁻ of
   theory.
4. The occupancy at that point is the measured occupancy, and the f″ value
   supports or contradicts the metal's identity.

This is the procedure described in the supplements of:

- Adv. Mater. 2023 (mmDNA 3D library), "Disambiguation of Metal Sites and f″
  Refinement", after Liu & Hendrickson, Acta Cryst. D 69, 1314 (2013).
- JACS 2023 (heterobimetallic pH paper), "Heavy Atom Identification" and
  "Anomalous Groups Refinement (f″)".
- Small 2025 (transmetalation), sections 1.5 and 1.6.
- Sci. Adv. 2026 (Ag nanoclusters), f″/Q method and Table S4.

It replaces the manual loop (refine, edit occupancy in Coot, import, refine
again), in which ordinary AutoRefine rounds would refine away the held
occupancy and B and let f′ drift.

## Command line

```
nasolve autofdp MODEL.pdb DATA.mtz
                [--output DIR] [--site SITE]... [--wavelength A] [--cif FILE]...
                [--element SITE=ELEMENT]... [--candidates ELEMENT[,ELEMENT...]]
                [--theory ELEMENT=FDP]... [--tolerance E] [--max-iterations N]
                [--cycles N] [--b-start surrounding|model|VALUE] [--b-step A2]
                [--occupancy-floor Q] [--map-sigma S]
                [--no-confirm-b] [--final] [--plan]
                [--no-coot] [--coot PATH] [--phenix-root PATH]
```

For 7CU_Ag, from the dataset folder:

```
nasolve autofdp AutoMR/run_001/AutoRefine/round_002/refined_001.pdb staraniso_alldata-unique.mtz
```

| Option | Default | Meaning |
|---|---|---|
| `MODEL.pdb` | required | the starting model, PDB format, for example an AutoRefine `refined_001.pdb`; never changed |
| `DATA.mtz` | required | reflections with anomalous data and Free-R flags. The two files may be given in either order |
| `--output` | `AutoFDP` beside the MTZ | folder in which a new `run_NNN` is made |
| `--site` | every metal site | metal sites whose B step 1 sets, e.g. `B:101`; repeatable. Every metal site is still refined in the f″ round |
| `--wavelength` | MTZ header | X-ray wavelength (Å), for an MTZ whose header has none; a disagreement with the header is recorded |
| `--cif` | none | restraint dictionary for a residue phenix.refine does not know; repeatable; also given to Coot |
| `--element` | element in the model | treat a site as another element, e.g. `B:101=HG` |
| `--candidates` | `AG,HG` plus the model's elements *(to confirm)* | elements compared for the identity check |
| `--theory` | cctbx Henke value at the wavelength | override theoretical f″, e.g. `HG=10.17`; repeatable. f′ stays the Henke value |
| `--tolerance` | 1.0 e⁻ | a site matches theory when \|refined f″ − theory\| ≤ tolerance |
| `--max-iterations` | 8 *(to confirm)* | f″ rounds before giving up |
| `--cycles` | 3 | phenix.refine macrocycles per f″ round |
| `--b-start` | `surrounding` | site B at the start: mean of the neighbouring nucleotides, the model's own (`model`), or a number |
| `--b-step` | 20 Å² *(to confirm)* | B change when the map says to adjust B instead of occupancy |
| `--occupancy-floor` | 0.05 *(to confirm)* | below this, the site is sent to identity review |
| `--map-sigma` | 3.0 *(to confirm)* | mFo-DFc peak or hole at the site that counts as significant |
| `--no-confirm-b` | confirm on | skip the B-confirmation round |
| `--final` | off | after convergence, refine XYZ, B, occupancy and f″ together (Adv. Mater. ending) |
| `--plan` | off | do the preflight and print what was read, without refining or writing anything |
| `--no-coot` | open Coot | do not open the f″ round's model and maps in Coot |
| `--coot` | configured Coot | Coot executable or installation |
| `--phenix-root` | configured Phenix | Phenix installation |

## What one run does

### 0. Preflight (implemented; fails before anything is written)

`plan_autofdp`. Everything AutoFDP needs comes from the two files:

- **Model**: an existing PDB-format file. Its metal sites (step 1) must be
  found, each with a known element.
- **MTZ** (`phenix.mtz.dump`):
  - the anomalous observations, F(+)/SIGF(+)/F(−)/SIGF(−) or
    I(+)/SIGI(+)/I(−)/SIGI(−). Without them f″ cannot be refined and AutoFDP
    stops;
  - the Free-R label (FreeR_flag, R-free-flags, FREE or FreeR). Flags are
    never generated. The test value is 1 for Phenix `R-free-flags` and 0 for
    the CCP4 labels, as AutoRefine uses for STARANISO data;
  - the wavelength from the MTZ header (the dataset that holds the data, as in
    AutoSol and AutoRefine), or `--wavelength` when the header has none.
- **Theory**: f′ and f″ of every metal element at that wavelength (cctbx Henke
  via the `phenix.python` beside `phenix.refine`), with `--theory` replacing an
  element's f″. Each theoretical f″ must be positive.
- `--cif` files must exist.

`--plan` stops here and prints the labels, Free-R flags, wavelength, theory,
sites and planned B changes. It runs `phenix.mtz.dump` and `phenix.python`
(read-only) but not `phenix.refine`, and writes nothing.

### 1. B start (implemented)

`prepare_start_model`:

1. **Source.** The model given on the command line. It is never changed.
2. **Copy.** The model is copied unchanged to `run_NNN/source_model.pdb`, and
   the copy's checksum is checked.
3. **Find the metals.** As in PostMR: NARestraints `find_heavy_atoms`, HETATM
   records with Z >= 21 outside nucleotides, without Li, Be, Na, Mg, Al, K or
   Ca. `--site` restricts the list. No metal, or an unknown `--site`, stops
   with an error. A metal with an ANISOU record is refused, because AutoFDP
   sets isotropic B only.
4. **Neighbours.** Found with the NARestraints metal-pair search
   (`restraints.metal_pairs`), in this order:
   1. the two nucleotides of the base pair the ion bridges
      (`find_metal_mediated_pairs`: base N, O, S or Se within 3.6 A of the ion
      in two nucleotides that are not chain neighbours, closest pair kept);
   2. otherwise every nucleotide with a base N, O, S or Se within 3.6 A;
   3. otherwise every nucleotide with any atom within 6 A.

   No nucleotide within 6 A is an error: remove the ion, or leave it out with
   `--site`.
5. **Set B.** The site's B becomes the mean B over all atoms of those
   nucleotides (first model), rounded to 0.01 (`--b-start surrounding`).
   `--b-start model` keeps the model's B; a number sets that B. Only the
   B-factor columns (61-66) of the metal lines change. Occupancy, coordinates
   and every other atom are unchanged. The result is `run_NNN/start_model.pdb`.
6. **Record.** `report.json` and `autofdp.log` hold, for each site:
   occupancy, B before and after, the neighbour rule used, and each
   neighbour's residue, mean B, atom count and closest contact.

`--plan` does steps 3-5 in memory and prints them, writing nothing. Each run
gets a new `run_NNN` folder; nothing is overwritten.

For 7CU_Ag `round_002/refined_001.pdb`, the Ag B:101 (B 379.99) bridges B:4 DU (residue mean
B 293.1) and A:12 DC (311.6).

AutoFDP expects to start from a model that has already been refined. It does
not do the initial coordinate refinement itself, and it does not place or
delete sites. Spurious or duplicate sites are removed by the user before
AutoFDP, as in the papers.

### 2. f″ round (implemented)

`fdp_round`. One phenix.refine run on `start_model.pdb`, in
`run_NNN/round_001/`:

- strategy `group_anomalous` only. Coordinates, real-space coordinates, B and
  occupancy are not refined, for any atom, and the X-ray/geometry and ADP
  weight optimisations are off;
- one anomalous group per metal site of the model (all of them, including any
  left out of step 1 with `--site`), each with f′ fixed at its theoretical
  value for the wavelength and f″ refined from 0;
- `--cycles` macrocycles (default 3), the MTZ's anomalous observations and
  Free-R flags, no experimental phases, and any `--cif` dictionaries.

The parameters (`autofdp.params`) and the command come from AutoRefine's
writers, so labels, Free-R handling and the Phenix 2.x data-manager format are
the same as in AutoRefine. phenix.refine checks them with `--dry_run` first.
Geometry restraints are not needed for an f″-only round, so no NARestraints or
secondary-structure file is passed; phenix.refine still needs a dictionary for
every residue it does not know, which `--cif` supplies (for example the PostMR
`ReadySet/prepared_model.ligands.cif`).

After the round, for each group:

- refined f″ (from the phenix.refine log), occupancy and B;
- deviation = refined f″ − theoretical f″, and whether it is within
  `--tolerance`;
- a check that coordinates, occupancy and B did not change (largest change
  recorded; a warning above 0.01 Å, 0.01 or 0.1 Å²).

The round's model and map coefficients open in Coot (2mFo-DFc and mFo-DFc with
`--auto`, and ANOM/PHANOM at 3σ when present), unless `--no-coot`. Coot's
files go to `run_NNN/CootGUI/`. A Coot problem never fails the run; the model
and map paths are printed instead.

Status `AUTOFDP_FDP_ROUND_DONE` when phenix.refine reported a refined f″ for
every group, otherwise `AUTOFDP_FAILED`.

Not yet implemented in this step: the mFo-DFc and anomalous map values at each
site (`map_value_at_site`).

### 3. Decision for each site

Let *q* be the occupancy, *f″ref* the refined value and *f″th* the theory.

- **Converged**: \|f″ref − f″th\| ≤ tolerance. The site is held from now on.
- **Otherwise**, the proposed occupancy is *q′ = q × f″ref / f″th* *(to
  confirm: the papers adjust occupancy by hand, in the direction given by the
  sign of the difference)*. Then:
  - Raising occupancy (*q′ > q*) when mFo-DFc at the site is already a
    significant hole (≤ −`--map-sigma`): lower B by `--b-step` instead.
  - Lowering occupancy (*q′ < q*) when mFo-DFc is already a significant peak
    (≥ +`--map-sigma`): raise B by `--b-step` instead.
  - *q′ > 1.0*: occupancy is set to 1.0 and the site is flagged "f″ above
    theory at full occupancy". Next round, B is lowered, if the map allows it.
    If f″ is still too high, the site goes to identity review (a heavier-f″
    element may fit better).
  - *q′ < `--occupancy-floor`*: identity review (a lighter-f″ element may fit
    better). For example, in the JACS paper an "Hg" site that stayed short of
    10.17 below occupancy 0.10 turned out to be Ag.
  - f″ that barely changes (by less than 0.2 × tolerance *(to confirm)*) after
    an occupancy change: flagged "unresponsive". The JACS paper took this as a
    sign of two metals or a wrong model.
- AutoFDP writes a model with the new occupancy/B values and runs the next f″
  round (`round_002`, ...) in the same run folder, for the sites that are not
  yet converged.
- After `--max-iterations` rounds, sites that have not converged go to review.

Sites are updated together in each round. Sites closer than 3.5 Å to each other
(for example the N3 and O4 metals of a bimetallic pair) have correlated f″
values, and they are reported as a coupled pair.

### 4. Identity check

After the first f″ round, for each candidate element *E* the implied occupancy
is *q × f″ref / f″th(E)*. An element is plausible if that occupancy is above
the floor and at most 1.0.

- If the model's element is plausible and no other element is: no message.
- If several elements are plausible (low-occupancy Hg looks like high-occupancy
  Ag): reported as ambiguous. The model's element is kept, and the papers'
  other evidence is listed as the user's decision: Fo-Fc after removing the
  site, Z (Hg 80 vs Ag 47), and the series trend.
- If the model's element is not plausible but another is: the run stops with
  review and suggests `--element SITE=E`.

AutoFDP never changes an element on its own. A new run with `--element`
renames the atom and residue in its copy of the model.

### 5. B confirmation (default on)

Following the Adv. Mater. supplement: one round with occupancy, f′ and f″
fixed at the converged values, refining only the B of the sites. If a site's B moves by more than 25% *(to confirm)*, the
site goes to review.

### 6. Final round (`--final`, default off)

The Adv. Mater. ending: XYZ, B, occupancy and f″ refined together, with f′
still fixed at theory. The JACS and Small supplements
stop at step 5. AutoFDP reports how far occupancy and f″ move in this round.

## Outputs

Each run is a new folder, `AutoFDP/run_NNN/` beside the MTZ by default (at the
same level as `AutoMR/` and `AutoSolMetal/` in a NASolve dataset), or inside
`--output`:

- `source_model.pdb` (unchanged copy of the given model) and
  `start_model.pdb` (metal B set), from step 1;
- `round_001/`: `autofdp.params`, `command.json`, `phenix.refine.preflight.log`,
  `phenix.refine.log`, and phenix.refine's outputs (`refined_001.pdb`,
  `refined_001_map_coeffs.mtz`, ...);
- `report.json`, with `workflow: "autofdp"` and schema 2:
  - `inputs`: the model, the MTZ and any `--cif` files, each as a portable
    reference (anchored to the dataset folder when the run is in
    `DATASET/AutoFDP/`), with checksum and size; the absolute path at creation
    is kept as provenance only;
  - `reflections`: labels, Free-R label and test value, wavelength and its
    source, resolution;
  - `settings`, `phenix_version`, and `theory` (per element: Henke f′, the f″
    used and its source, Henke f″);
  - `b_start`: the two model files and, per site, occupancy, B before and
    after, and the neighbours used;
  - `rounds`: per round, the folder, strategy, starting groups, cycles, R-work
    and R-free, the held-parameter check, the refined model and map
    coefficients, and per site the occupancy, B, fixed f′, refined f″, theory,
    deviation and whether it is within tolerance;
  - `status`, `message` and `steps`.

  Planned, in addition: mFo-DFc and anomalous map σ at each site, the decision
  per round, the identity check, and the recommended model.
- `autofdp.log`: the same in human-readable form;
- `CootGUI/`: Coot's log, backups and map script;
- planned: `fdp_table.tsv`, one row per site per round, for plotting.

PDB files do not store f″; the refined values are in `report.json`.

## Status and exit codes

| Status | Exit | When |
|---|---|---|
| `AUTOFDP_PREPARED` | - | in `report.json` between step 1 and the end of the round |
| `AUTOFDP_FDP_ROUND_DONE` | 0 | step 2 done: one f″ round with a refined f″ for every group (until occupancy adjustment exists) |
| `AUTOFDP_CONVERGED` | 0 | every site within tolerance and B confirmation passed |
| `AUTOFDP_REVIEW` | 2 | any site not converged, hit an occupancy bound, unresponsive, coupled, ambiguous in identity, or with B drift; per-site reasons are in the report |
| `AUTOFDP_FAILED` | 1 | phenix.refine failed or reported no f″; the run folder is kept |
| preflight error | 1 | before any folder is written |

## Interaction with other commands

- AutoFDP is independent of AutoMR, PostMR and AutoRefine: it reads no run
  report or checkpoint and writes nothing outside its own run folder. A
  NASolve dataset's MTZ and an AutoRefine `refined_NNN.pdb` are the usual
  inputs, but any model and MTZ with anomalous data work.
- Refine Doctor is unchanged.
- Campaigns do not run AutoFDP.

## AutoFDP never

- regenerates Free-R flags;
- refines coordinates, or refines f′, during f″ rounds;
- changes the given model or MTZ, or atoms other than the sites in its copy;
- changes an element without `--element`;
- places or deletes sites;
- writes outside its own run folder, or overwrites an earlier run.

## Questions for the lab before implementation

1. Is occupancy adjusted proportionally (*q × f″ref / f″th*) or by hand in
   steps?
2. Is f″ restarted from 0 in every round, or from the previous value?
3. "B roughly the same as its neighbours" is implemented as the atom-weighted
   mean B of the bridged pair (with the fallbacks above). Should it be all
   atoms of those nucleotides (current), only their base atoms, or also the
   stacked neighbours?
4. When should B change instead of occupancy? Is the threshold `--map-sigma`
   3σ, and is `--b-step` 20 Å² a fixed step, or should it scale with
   resolution?
5. Should the run end at ±1 e⁻ (JACS, Small) or with the all-parameter
   refinement (Adv. Mater.)?
6. Which theory values should be used? cctbx Henke gives Ag 2.03 e⁻ at
   1.00743 Å; the papers quote Ag 2.06 and Hg 10.17 at 12300 eV. Tabulated
   values ignore white-line effects at an absorption edge (the JACS supplement
   measured Au f″ 11 against about 10 at the Au L-III edge), hence `--theory`.
7. Are the occupancy floor (0.05), the iteration limit (8) and the B-drift
   limit (25%) reasonable?
8. Should coupled sites (bimetallic pairs) be calibrated together, or one at a
   time with the other held?
