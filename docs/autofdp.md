# AutoFDP: occupancy from f″ refinement

**Status: step 1 implemented, later steps specified only (2026-11-11).**
`nasolve autofdp` runs step 1 (B start, below): it copies an AutoRefine result
and sets each metal's B to the mean B of its neighbouring nucleotides, as a
new checkpoint, then stops. The f'' rounds (steps 2-6) are documented stubs in
`src/nasolve/autofdp.py` that raise `NotImplementedError`; options that belong
to them are accepted but have no effect yet, and `--element` is refused. The
defaults marked *(to confirm)* follow my reading of the papers below and have
not yet been checked by the lab.

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
nasolve autofdp [RUN] [--from CHECKPOINT]
                [--site SITE]... [--element SITE=ELEMENT]...
                [--candidates ELEMENT[,ELEMENT...]] [--theory ELEMENT=FDP]...
                [--tolerance E] [--max-iterations N] [--cycles N]
                [--b-start surrounding|model|VALUE] [--b-step A2]
                [--occupancy-floor Q] [--map-sigma S]
                [--no-confirm-b] [--final] [--plan]
```

| Option | Default | Meaning |
|---|---|---|
| `RUN` | active workspace run | AutoMR run folder with a completed PostMR |
| `--from` | current checkpoint | starting model: an AutoRefine checkpoint ID or bookmark |
| `--site` | every metal site | restrict to these metal sites, e.g. `B:101`; repeatable |
| `--element` | element in the model | treat a site as another element, e.g. `B:101=HG`; renames the atom in a new checkpoint |
| `--candidates` | `AG,HG` plus the model's elements *(to confirm)* | elements compared for the identity check |
| `--theory` | cctbx Henke value at the wavelength | override theoretical f″, e.g. `HG=10.17` |
| `--tolerance` | 1.0 e⁻ | convergence: \|refined f″ − theory\| ≤ tolerance |
| `--max-iterations` | 8 *(to confirm)* | f″ rounds before giving up |
| `--cycles` | 3 | phenix.refine macrocycles per f″ round |
| `--b-start` | `surrounding` | site B at the start: mean of the neighbouring nucleotides, the model's own (`model`), or a number |
| `--b-step` | 20 Å² *(to confirm)* | B change when the map says to adjust B instead of occupancy |
| `--occupancy-floor` | 0.05 *(to confirm)* | below this, the site is sent to identity review |
| `--map-sigma` | 3.0 *(to confirm)* | mFo-DFc peak or hole at the site that counts as significant |
| `--no-confirm-b` | confirm on | skip the B-confirmation round |
| `--final` | off | after convergence, refine XYZ, B, occupancy and f″ together (Adv. Mater. ending) |
| `--plan` | off | print the sites and B changes without writing anything (later: the full preflight) |

## What one run does

### 0. Preflight (fails before anything is written)

- The run has a completed PostMR. The starting checkpoint exists and is usable.
- The observations have an anomalous quartet (F(+), SIGF(+), F(−), SIGF(−), or
  the intensity equivalent) and a Free-R label. Flags are never regenerated.
- The wavelength comes from the MTZ header, with the AutoSol record as a
  fallback (as in AutoRefine). Without a wavelength, stop.
- The sites are the PostMR metal sites plus nucleotide heavy-atom candidates,
  or those named with `--site`. Each must occur exactly once in the starting
  model, with a known element.
- Theoretical f′ and f″ are calculated for every site element and candidate
  element at the wavelength (cctbx Henke via `phenix.python`), unless
  `--theory` overrides them. Each theoretical f″ must be positive.
- The Phenix version is recorded, and `phenix.refine --dry_run` must accept the
  f″-only parameters.

`--plan` stops here and prints the sites, elements, theoretical values,
starting occupancy and B, and the settings.

### 1. B start (implemented)

`nasolve autofdp RUN --from refine-NNN`:

1. **Source.** The starting checkpoint (`--from`, default the current one)
   must be an AutoRefine refinement (`kind: refinement`) that is usable. Its
   model checksum is verified. A PostMR or imported checkpoint is refused with
   the list of AutoRefine checkpoints. This matters because rounds that end in
   review are not selected as current, so the current checkpoint is often
   still `postmr`.
2. **Copy.** The model is copied unchanged to
   `AutoFDP/fdp_NNN/source_model.pdb`, and the copy's checksum is checked.
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
   and every other atom are unchanged. The result is
   `AutoFDP/fdp_NNN/start_model.pdb`.
6. **Checkpoint.** The start model is registered as checkpoint `autofdp-NNN`
   under the source checkpoint. It is a model edit outside phenix.refine, like
   an imported model, so it has kind `manual`, with recipe `AutoFDP/b-start`,
   status REVIEW, the parent's observations, phases and restraints, and an
   `autofdp` block listing the edits. `nasolve show RUN --checkpoint
   autofdp-NNN` shows it with its parent's maps. The current checkpoint is not
   changed.
7. **Record.** `AutoFDP/fdp_NNN/report.json` (status `AUTOFDP_PREPARED`) and
   `autofdp.log` hold, for each site: occupancy, B before and after, the
   neighbour rule used, and each neighbour's residue, mean B, atom count and
   closest contact.

`--plan` does steps 1-5 in memory and prints them, writing nothing. Each run
gets a new `fdp_NNN` folder and checkpoint; nothing is overwritten.

For 7CU_Ag `refine-002`, the Ag B:101 (B 379.99) bridges B:4 DU (residue mean
B 293.1) and A:12 DC (311.6).

AutoFDP expects to start from a model that has already been refined. It does
not do the initial coordinate refinement itself, and it does not place or
delete sites. Spurious or duplicate sites are removed by the user before
AutoFDP, as in the papers.

### 2. f″ round

One phenix.refine run, recorded as a normal AutoRefine round
(`AutoRefine/round_NNN`, checkpoint `refine-NNN`, recipe `AutoFDP/fdp`):

- strategy `group_anomalous` only. No coordinates, B or occupancy are refined,
  for any atom;
- one anomalous group per site, with f′ fixed at its theoretical value and f″
  refined from a starting value of 0;
- `--cycles` macrocycles, the same observations, Free-R flags, restraints and
  phases as the starting checkpoint.

For each site, the round records:

- refined f″;
- the mFo-DFc value at the site, in σ;
- the anomalous difference map value at the site, in σ (reported only).

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
- AutoFDP writes a checkpoint with the new occupancy/B values (a model-edit
  checkpoint like step 1's) and runs the next f″ round, for the sites that are not yet
  converged.
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
renames the atom and residue in a new checkpoint.

### 5. B confirmation (default on)

Following the Adv. Mater. supplement: one round with occupancy, f′ and f″
fixed at the converged values, refining only the B of the sites (recipe
`AutoFDP/confirm-b`). If a site's B moves by more than 25% *(to confirm)*, the
site goes to review.

### 6. Final round (`--final`, default off)

The Adv. Mater. ending: XYZ, B, occupancy and f″ refined together, with f′
still fixed at theory (recipe `AutoFDP/final`). The JACS and Small supplements
stop at step 5. AutoFDP reports how far occupancy and f″ move in this round.

## Outputs

Every phenix.refine run is an ordinary AutoRefine round and checkpoint, so
`nasolve checkpoints list` and `nasolve show RUN --checkpoint ID` work on every
step. The AutoFDP record goes to `AutoMR/run_NNN/AutoFDP/fdp_NNN/`, beside
`AutoRefine/` and `RefineDoctor/`:

- `source_model.pdb` (unchanged copy of the AutoRefine result) and
  `start_model.pdb` (metal B set), from step 1;
- `report.json`, with `workflow: "autofdp"` and schema 1, containing:
  - the settings, the Phenix version, and the wavelength with its source;
  - the theory table: f′ and f″ per element, and whether each value came from
    cctbx Henke or from `--theory`;
  - the sites and their iterations. For each site and round: the checkpoint,
    element, occupancy, B, fixed f′, f″ start, refined f″, theory, difference,
    mFo-DFc and anomalous map σ at the site, and the decision;
  - the identity check, the final status with a reason per site, and the
    recommended checkpoint.
- `autofdp.log`: the human-readable progress.
- `fdp_table.tsv`: one row per site per round, for plotting.

PDB files do not store f″. Each AutoFDP checkpoint therefore records its
anomalous groups (selection, f′, f″, and which values were refined) in its
checkpoint metadata.

## Status and exit codes

| Status | Exit | When |
|---|---|---|
| `AUTOFDP_PREPARED` | 0 | step 1 done (the only step implemented so far) |
| `AUTOFDP_CONVERGED` | 0 | every site within tolerance and B confirmation passed |
| `AUTOFDP_REVIEW` | 2 | any site not converged, hit an occupancy bound, unresponsive, coupled, ambiguous in identity, or with B drift; per-site reasons are in the report |
| `AUTOFDP_FAILED` | 1 | phenix.refine failed; completed rounds and checkpoints are kept |
| preflight error | 1 | before any folder or checkpoint is written |

## Checkpoints

AutoFDP never changes the current checkpoint, the same as Refine Doctor. On
`AUTOFDP_CONVERGED`, it bookmarks the final checkpoint as `autofdp-NNN` and
prints `nasolve checkpoints use RUN autofdp-NNN`.

## Interaction with other commands

- `nasolve autorefine --from <AutoFDP checkpoint>` should keep the calibrated
  f′/f″ fixed instead of refining both again from 0. This needs a small
  AutoRefine change when AutoFDP is implemented.
- Refine Doctor is unchanged.
- Campaigns do not run AutoFDP.

## AutoFDP never

- regenerates Free-R flags;
- refines coordinates, or refines f′, during f″ rounds;
- changes atoms other than the sites;
- changes an element without `--element`;
- places or deletes sites;
- overwrites or deletes earlier rounds or checkpoints;
- changes the current checkpoint.

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
