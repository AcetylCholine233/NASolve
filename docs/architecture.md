# NASolve architecture

NASolve orchestrates external crystallographic programs; it does not reimplement them.

## Runtime precedence

Phenix is resolved in this order:

1. `--phenix-root`
2. `NASOLVE_PHENIX_ROOT`
3. saved user configuration
4. current `PATH`
5. validated installations in standard platform locations
6. a clear request to run `nasolve configure phenix PATH`

When multiple standard installations exist, NASolve selects the highest validated version.
An automatically discovered path is saved and revalidated on use. Each structure-solution report will record the exact
Phenix version and executable paths used.

Coot uses the same pattern independently: a one-run path, `NASOLVE_COOT`,
saved configuration, `PATH`, then bounded standard locations. Validation checks
both `coot --version` and the embedded Python API required by headless mutation.

## AutoMR input contract

`nasolve automr DATASET` uses one shared `nasolve.txt` schema for standard and
nonstandard models. The file records user intent; computed paths, checksums,
model statistics, and resolved ligand codes belong in the run log and JSON
report.

Standard example:

```ini
[automr]
mode = standard
frame = W
pair = D:T
```

Nonstandard example:

```ini
[automr]
mode = nonstandard
model = custom_model.pdb
sequence_file = construct.fasta
mirror = false

[mutations]
A:3 = 5IU
```

`sequence_file` accepts chain-labelled FASTA or `CHAIN = SEQUENCE` records.
Inline `[sequences]` remains an exclusive alternative. The effective sequence
is frozen inline in the run snapshot, with the source path and checksum retained
in the report.

If `nasolve.txt` is absent, the command infers a nonstandard run from exactly
one top-level PDB, or takes a standard frame from `-W`, `-3GBI`, or `--frame`.
It then writes the minimal canonical file. Existing input is never overwritten.
CLI frame and pair values override file values and are captured in the frozen
run snapshot.

Current standard models are found in repository `MR_frames/`, through
`NASOLVE_MR_FRAMES`, or with `--frames-dir`. Approved models will later move
into package data for installed releases. `W` and `5W6W` resolve to the same
frame; absence of an explicit standard frame never defaults to W.

The model catalogue is organized by frame and ordered pair:

```text
MR_frames/
├── 5W6W/
│   ├── A_G.pdb
│   ├── C_G.pdb
│   └── D_F.pdb
└── 3GBI/
    └── C_C.pdb
```

Names use `FIRST_SECOND.pdb`, where each component is either a configured alias
or a literal NARestraints ligand code. Selection compares the resolved ligand
codes, not just filename spelling. Thus `D_F.pdb` matches either `D:F` or
`1AP:DF`. The deterministic fallbacks are `5W6W/C_G.pdb` and
`3GBI/C_C.pdb`. An exact pair model needs no pair mutation; a fallback records
the requested post-MR mutation. Missing fallbacks, malformed names, and
duplicate exact matches stop the run. An unrelated catalogue entry with an
unrecognized ligand token is skipped and recorded as a warning; it does not
block a valid requested pair. Requesting that unknown token still fails normal
ligand validation.

Every preflight receives an isolated directory:

```text
DATASET/AutoMR/run_001/
├── nasolve.input.txt
├── automr.log
├── report.json
└── Model/
    ├── source_model.pdb  # only when mirroring
    ├── input_model.pdb
    └── assessment.json
```

Adding `--execute` continues that same frozen run with an isolated Phaser
directory:

```text
DATASET/AutoMR/run_001/Phaser/
├── phaser.eff
├── phaser.log
├── PHASER.1.pdb
├── PHASER.1.mtz
├── mr_solution.pdb
└── mr_solution.mtz
```

The generated PHIL defines exactly one explicit ensemble named
`nasolve_model`; NASolve does not also use Phaser's `model` shortcut. It sets
`use_hetatm = True` and an explicit model RMSD of 1.0 A for compatibility with
Phenix builds that require either RMSD or fractional identity. It supplies one
NA composition using the model checker's polymer-residue count (`nres`), and
searches for one copy in H3/R3. The guarded
P1 shunt sets both composition and search copies to three. Phaser's original
outputs remain untouched, while stable `mr_solution` copies identify the
result selected for the next layer.

Without mirroring, the copied model is checksum-verified, so HETATM records and
all other bytes are preserved. With `--mirror` (or `mirror = true`), NASolve
copies the untouched source into the run, calls
`restraints.mirror_pdb.mirror_pdb`, and assesses/checksums the transformed
`input_model.pdb` that Phaser will consume. The mirror is rejected if atom,
polymer-residue, chain, or HETATM inventories change unexpectedly.

Residues written as HETATM whose residue names occur in the
NARestraints ligand library are counted as modified polymer residues, not as
detached ligands. This keeps chain lengths correct for modified nucleotides
while retaining separate counts for all HETATM atoms and nonpolymer hetero
residues. Full sequences must match the resulting model-chain residue count,
and explicit mutation targets must exist. Mutation precedence is full sequence,
standard pair, then explicit mutation.

### Dataset and symmetry discovery

Dataset discovery accepts the usual autoPROC/STARANISO naming variants rather
than renaming raw inputs. It prefers a unique MTZ whose punctuation-insensitive
name contains both `staraniso` and `alldata`, falling back to the only top-level
MTZ when there is exactly one. Metadata is the unique case-insensitive
`Data_1*.cif`; multiple candidates stop the run.

The symmetry gate applies only to standard W/3GBI recipes. NASolve compares:

1. the MTZ file symbol, number, and matrix symbol reported by
   `phenix.mtz.dump`;
2. every explicit space-group value in `Data_1*.cif`; and
3. only the final (bottom-most) `Spacegroup name` entry in `summary.html`.

Earlier indexing and POINTLESS alternatives in the HTML are ignored. `H 3`,
`H3`, `R 3`, `R3`, and the MTZ matrix form `R 3 :H` are accepted as one
equivalent class and use one MR copy. The exact MTZ label is retained in the
report as a required final-output patch after structure solution. A genuine P1
standard run stops by default; `--allow-p1-standard` records a strong red flag
and plans three MR copies. Other groups and disagreements between authoritative
sources stop the run. Nonstandard mode does not apply this symmetry rule.

## Pipeline gates

1. Validate and freeze AutoMR input; optionally mirror the search model; assess
   the exact model supplied to Phaser.
2. Run Phaser in an isolated `Phaser/` directory without stripping heteroatoms.
3. Classify the best TFZ: `>= 8.0` pass; `7.0–7.99` review; `< 7.0` fail.
4. `nasolve postmr RUN` applies supported site changes, builds the restraint
   stack, and runs ReadySet with hydrogens disabled.
5. If PostMR finds nucleotide-bound I, Br, or Se, `nasolve autosol RUN` performs
   guarded MR-SAD phasing and validates at least one heavy-atom site against
   the model under crystallographic symmetry.
6. `nasolve autorefine RUN` performs one quiet five-cycle default refinement,
   validates the output model and statistics, and appends an immutable
   checkpoint.
7. `nasolve refine-doctor RUN` audits a selected refinement and, when safe,
   runs a finite set of sibling diagnostic branches without changing the
   current checkpoint.

The default command still stops after step 1. `nasolve automr ... --execute`
runs step 2 and applies step 3. A passing run exits 0; TFZ 7.0–7.99 is retained
for review and exits 3; failure or missing output exits 4.

## PostMR contract

PostMR is deliberately a separate command so an accepted Phaser result can be
inspected and rerun deterministically before refinement. It refuses to
overwrite an existing `PostMR/` directory. `MR_SUCCESS` is accepted;
`MR_REVIEW` requires an explicit override; other statuses stop.

The W-frame manifest assigns requested pair roles to `A:12` and `B:4`.
Canonical DNA/RNA changes use a generated headless Coot script. Supported
modified nucleotides use a canonical-parent mutation followed by dictionary
monomer construction, overlap, and residue replacement. Coot runs with its
backup directory redirected beneath the run and removes hydrogens before
writing. Full sequences expand over the frozen per-chain residue inventory;
application precedence is sequence, standard pair, then explicit mutation.
Every ordinary base change restores the pre-Coot sugar/phosphate coordinates.
For mutations outside the curated registry, `ligand_definition` requires one
NARestraints record and derives a compatible Coot construction parent from
`Sugar Type` and `Base Analog`. `ligand_dictionary` requires a local
`ligands/CODE.cif`; no network retrieval occurs. The generic dictionary check
verifies `_chem_comp.id`, while curated entries keep their explicit topology
checks and take precedence. This permits OHU mutation with the included CIF
without claiming general topology or atom-mapping inference. The inferred Coot
parent is not necessarily the component's official CCD parent.

For mirrored runs, canonical targets are translated to the NARestraints L-side
codes (`0DA`, `0DC`, `0DG`, `0DT`, and RNA equivalents). Exact mirrored models
therefore are not accidentally mutated back to D chemistry. A required
mirror-side sequence change stops until the guarded
unmirror/Coot/remirror construction path is implemented.

For supported dictionary modifications, Coot also writes a canonical-parent
snapshot. NASolve copies the coordinates, occupancies, and B factors of every shared
atom back into the replaced residue and requires the complete sugar/phosphate
atom set to survive. Thus dictionary overlap positions only the new chemistry;
it cannot rotate or curl the existing phosphate. The restoration inventory is
recorded in the run report.

Dictionary resolution follows mutation actions. Curated components present in
the prepared model are also collected even when unchanged, so their CIFs and
`component_identity` entries survive a no-op mutation plan. This does not imply
automatic resolution of every pre-existing noncurated residue. Generic
definitions add no inferred sulfur substitutions or ring-substituent rules.

The registry also declares exact parent/target atom substitutions. `DE`, `DF`,
and `S6G` place sulfur along the canonical `C-O` vector, use the transformed
dictionary model only for the reviewed C-S bond length, and inherit occupancy
and B factor from the replaced oxygen. This prevents an otherwise valid Coot
graph overlap from placing the one new sulfur atom on the wrong side of a base.
`C38` and `5IU` use the same guarded construction with parents `DC` and `DT`.
Their iodine is placed along the external bisector defined by the mapped
`C4-C5-C6` ring atoms, using the C-I distance from the dictionary's ideal
coordinate set and the parent `C5` occupancy and B factor. The transformed Coot
monomer is not trusted for the halogen bond length.

For W/5W6W, NARestraints consumes the packaged `Std_padd.txt` specification
(`A 11:13` / `B 5:3`) and remains the sole source of stacking restraints. The
portable 5W6W secondary-structure template contains 17 non-overlapping
scaffold base pairs, no stacking pairs, and no historical input paths, cell,
map, output, or GUI settings.

`--modified-pairs-only` generates a model-derived restraint selection in either
standard or nonstandard mode. After all coordinate edits,
the NARestraints guesser runs with `allow_noncanonical=True`. NASolve defines
modified sites by noncanonical nucleotide residue identity and retains a
candidate when either partner is modified. It writes pair restraints with
stacking disabled. When a project also supplies a scaffold EFF, the generated
PHIL is authoritative: matching scaffold base-pair blocks are removed and the
remaining project overlay is retained. A zero-pair result is a successful
no-op; warnings, candidate counts, retained pairs, and removed overlaps remain
in the report.

Problematic Phenix components are explicit curated-library entries. User token
`E` maps to the laboratory PDB-compatible `DE`, not official CCD `8RO`, whose
topology describes a different sulfur-ring compound. `F` maps to `DF`, with
official deposition identity `A1AAZ`; `D` maps to official `1AP`; and `Q` maps
to official `S6G` with canonical parent `DG`. DNA `U` maps to official `DU`
(2'-deoxyuridine), whose packaged CIF is the unchanged wwPDB definition from
Phenix 2.2.1-6174 and whose canonical parent is `DT`. The registry
defines canonical parents and required/forbidden sulfur bonds, which are
validated before Coot or ReadySet runs. ReadySet receives the curated CIF,
executes with `actions.hydrogens=False`, and is run with `cwd` set to its own
directory because Phenix 1.20.1 accepts but ignores `input.output_dir`. The
output is rejected if it contains hydrogen coordinates. A phosphate check
removes verified internal OP3/O3P additions before comparing total and HETATM
atom counts to the prepared input; other count changes remain errors.

`phosphate.sanitize_phosphates` also runs before NARestraints. It scopes checks
to nucleotide-like residues with OP3/O3P and to sites repaired in the first
pass. An incoming O3'-P contact must be unique within a model/chain/TER segment,
with a 1.2–2.1 A separation. Remaining P, OP1/O1P, OP2/O2P and O5'/O5* atoms
must be unambiguous. Their P-O distances use the same broad bounds, and their
O-P-O angles must be between 60 and 160 degrees. These are coarse corruption
guards, not chemical target values. Reference-known incoming and outgoing
links must survive. OP3/O3P requires explicit per-site user consent, including
at termini. Selecting an annotated recipe is also explicit consent for exactly
its declared `chemistry.terminal_phosphate_sites`: the W/5W6W card declares
D:1, confirmed by Simon as part of the design on 2026-09-15. This is resolved
before AutoMR freezing, not inferred later from coordinates or residue numbers.
An explicit dataset `allow_op3_sites` replaces the entire recipe list, including
an empty override. Unrequested unlinked OP3 stops for review. A request must not
permit an extra OP3 on an internal phosphate. Ambiguous explicit or symmetry
links require inspection. Only verified extra atoms and associated
ANISOU/SIGATM/SIGUIJ/CONECT records are removed; surviving coordinates,
occupancies, and B factors are untouched.

Atom serials need not be unique for this coordinate cleanup: selected atoms
are identified by their validated coordinate records, while ANISOU/SIGATM/
SIGUIJ records use exact atom name, alternate, residue identity, model and
TER segment. Unrelated records sharing the same serial are preserved. Because
CONECT identifies atoms only by serial, affected bonds must have unique
endpoints before they can be interpreted or edited. Cleanup does not renumber
the raw or surviving coordinate records.

PostMR's additive `phosphate_cleanup.before_restraints` and `after_readyset`
reports record checked sites, removed/retained atoms, connectivity, and measured
geometry. Raw ReadySet output remains at `readyset.updated_model`; the new
`readyset.phosphate_checked_model` points to the separately checked copy used
for the final model. Immutable earlier runs require no migration and are never
rewritten.

The first parameterized monomer override uses selected-site `NASnoOP3`
modifications for reviewed 1AP chemistry. The original integration/validation
diary is preserved under `history/1ap-phosphate-integration.md`; it is no longer
the active contract.
New 1AP runs freeze authoritative refinement/view artifact lists with hashes,
so ReadySet cannot silently supersede the reviewed 1AP CIF or omit the
modification at checkpoint creation. Incoming/outgoing profile links and
explicit OP3 consent are checked after preparation, around refinement, before
Doctor triage and on manual import. Existing histories are never rewritten.

Coot/model editing, component CIFs, NARestraints, and ReadySet have separate
roles. Coot establishes coordinate identity and placement; CIFs establish
intra-residue topology; NARestraints establishes selected inter-residue
geometry. ReadySet is the uniform post-edit normalization and validation gate,
even where Phenix refinement could technically consume the edited PDB and CIF
directly.

NARestraints v1.1.1 has a legacy `DF` mapping of canonical `O2` to `S2`; the
A1AAZ-derived component uses `S1`. The PostMR adapter patches a copied residue
record in process memory, restores the builder's original loader in a `finally`
block, and records the correction. The installed workbook is never mutated,
and a future upstream `S1` mapping requires no correction.
The same adapter trims the exact legacy `S6G C5` value `C5 ` to `C5`.

The ReadySet model is scanned by element rather than by a closed residue list.
Iodine, bromine, and selenium in nucleotide-like residues are the initial
anomalous trigger set. The PDB element column is authoritative; a conservative
atom-name fallback is recorded when the column is blank. This allows a new
modified residue to trigger the later AutoSol branch without first adding its
residue code to an AutoSol allowlist.

## AutoSol contract

AutoSol is conditional and deliberately separate from PostMR. It requires a
`POSTMR_READY` report containing one supported heavy-atom element and refuses
to overwrite an existing `AutoSol/` directory. The current trigger elements
are iodine, bromine, and selenium; the element registry is intentionally open
to reviewed additions.

When a standard-frame catalogue supplies `seq_base.txt`, AutoMR freezes a
checksummed copy under `Model/` and AutoSol copies that run-local snapshot into
its input. Legacy runs recover the exact declared frame directory without
requiring the old pair-model filename to remain present. For a nonstandard
model, complete sequences recorded in the frozen run may be written to the same
input instead. The wavelength is read from the MTZ header (the dataset that
holds reflection data); a present `summary.html` must agree within 1e-4 A.
`phenix.mtz.dump` must confirm a complete anomalous intensity or amplitude
quartet before execution.

The AutoSol request is MR-SAD guided by the original Phaser solution PDB, not
the modified PostMR model. This avoids using the modeled heavy atom to validate
its own anomalous signal. NASolve supplies the heavy-atom element and
wavelength, leaves the number of sites unspecified, uses every logical CPU,
and explicitly sets both `build=False` and
`phase_improve_and_build=False`. After execution it reads `autosol.eff` and
rejects the run unless both build controls remain false, every `sites` value is
unset, and the requested processor count is effective.

Output discovery is content-based. The refinement input must contain a full
Hendrickson-Lattman phase quartet and anomalous observations. The AutoSol HA
PDB is compared with every expected nucleotide heavy atom of the same element.
For H3/R3, matching expands the three point operators, R-centering
translations, and neighboring unit cells; P1 uses lattice translations. The
branch automatically accepts phases when at least one expected atom has a
same-element HA site within 4 A. A nearest site between 4 and 8 A produces an
`AUTOSOL_REVIEW` result and requests inspection. A more distant or absent site,
or an AutoSol execution/output failure, produces `AUTOSOL_WARNING`. Review and
warning outcomes do not block the ordinary MR refinement path, but their phase
files are not automatically approved for use. All sites and match attempts are
retained in the report, while extra unmatched sites do not block this phase of
the pipeline.

AutoSol also records an optional `outputs.density_modified_map` portable
reference and SHA-256, plus `density_modified_map_status` and a diagnostic.
Discovery prefers a unique `overall_best_denmod_map_coeffs.mtz`, then a unique
`resolve_histograms*.mtz`, excluding TEMP/PDS products. Missing, ambiguous, or
unreadable view maps do not change phase acceptance. This viewing product is
distinct from the HL-bearing `refinement_data` input.

## Metal-search AutoSol contract

`nasolve autosol RUN --metal ELEMENT` is a separate, optional MR-SAD search for
metal ions, run directly after Phaser. It requires an AutoMR run with an
accepted Phaser result (best TFZ >= 8; 7 <= TFZ < 8 only with
`--allow-mr-review`) and may also run after PostMR, because it always uses the
original Phaser solution PDB and never a modified model. Inputs follow the
nucleotide branch: checksum-verified STARANISO reflections, the MTZ-header
wavelength, and a complete anomalous quartet confirmed by `phenix.mtz.dump`.
The sequence is, in order, `--sequence FILE`, the frozen frame or run sequence
(a missing or altered frozen sequence is still an error), or the nucleotide
chains of the Phaser model (residues with C1'; curated modified residues count
as their parent, DU as T, others are left out with a warning). The 3GBI
catalogue has no `seq_base.txt`, so 3GBI runs use the Phaser model. The report
records `sequence_source` and `sequence_warnings`. All checks complete before
any directory is created, so a refused search leaves nothing behind.

The first named metal is `autosol.atom_type`; additional metals go to
`autosol.phasing.mad_ha_add_list`. `autosol.phasing.partpdb_rms` defaults to
1.0 A (`--rmsd`). Sites are left unspecified and both
`model_building.build` and `phase_improve_and_build` are false. After a
successful exit, `autosol.eff` must confirm those settings and the requested
`partpdb_rms`; a mismatch fails closed.

Each search allocates a new `DATASET/AutoSolMetal/run_NNN`, beside `AutoMR`
and the input data. Outcomes are
`METAL_SITES_FOUND` (at least one site of a requested element in a unique
`overall_best_ha_pdb.pdb`), `METAL_NOT_FOUND` (no such site; a normal result),
and `METAL_AUTOSOL_FAILED` (non-zero exit, ambiguous outputs, or unapplied
settings). For each found site the search records the closest atom of the
three nearest Phaser-model residues, using H3/R3/P1 symmetry when configured
and same-asymmetric-unit distances otherwise.

The full record is the search directory's own `report.json`
(`workflow = "metal_autosol"`, `schema_version` 2, `source_run`) with SHA-256
checksums and references anchored to the search directory (its log, settings
file, sequence and HA model) or the dataset (reflections, Phaser model and
map), so it validates after the dataset moves. The AutoMR run's
`report.json["metal_autosol"]` lists the searches made from it (`searches`,
`latest`, each with a checksummed reference to the search report). Schema-1
`attempts` written inside the AutoMR run by development builds remain readable
by the metal view. The run's `stage`,
`status` and `message` are unchanged and `model_changed` is false: sites are
for inspection only and never enter the model or AutoRefine automatically.
The `metal` Coot view opens the Phaser model and the HA model over the Phaser
MTZ in `CootGUI/metal/` of the directory that was opened (the search directory
after `nasolve autosol --metal`), refusing references that fail their
checksums. Its startup script also adds a NASolve toolbar action (and
`nasolve_merge_metal_sites()`) that merges the open HA molecule into a copy of
the Phaser model, caps occupancies at 1.00, and writes a new
`<model>-merged[-N].pdb` in the source AutoMR run's `Phaser/` directory without
overwriting, then appends it to that run's `report.json["metal_merge"]` (run-
anchored, SHA-256) and selects it. PostMR starts from the selected merged model
unless `--original-model` is given, refuses one that is missing or changed, and
records `input_model_source`; AutoSol metal searches keep using the original
Phaser model. Live validation against Phenix 2.2.1-6174 AutoSol is still pending;
tests use a simulated AutoSol that accepts only the full PHIL paths above.

## AutoRefine contract

AutoRefine separates immutable scientific inputs from evolving coordinates.
Each round inherits an improved model from one selected parent checkpoint, but
always retains the original STARANISO observations and original Free-R array.
Refinement-generated MTZ files are outputs and evidence; they cannot become the
observation source of a child by accident.

The root checkpoint is the ReadySet model. Its restraint set prefers
ReadySet's generated combined ligand CIF over the individual curated CIFs that
were supplied to ReadySet, while retaining the NARestraints PHIL and project
EFF files. An approved AutoSol file is a separate phase source. Explicit
reflection file names and labels prevent duplicate observations in the AutoSol
file from competing with the authoritative STARANISO array.
AutoRefine derives one fail-closed reflection-selector policy from the
discovered Phenix version before allocating a round. Phenix 1.20.x uses
`legacy-explicit`: the command retains exact observation, Free-R, and optional
HL phase file/label selectors, but the unsupported Data Manager block is
omitted. Phenix 2.1.x and 2.2.x use `data-manager-file-scoped`: the same explicit
selectors are retained and the parameter file additionally binds each label
set to its MTZ through the Data Manager. Unknown, malformed, and unvalidated
families are rejected rather than mapped by version ordering. The exact
version and selected mode are frozen in the round report and checkpoint, and
Refine Doctor passes the same policy to every trial.
AutoRefine hashes that phase source once immediately before use, captures its
filesystem identity, and checks the identity again after all external work and
output hashing. Replacement, modification, or deletion is retained as a failed
checkpoint, never selected as current, so the next attempt receives a new
immutable round number.

The default recipe is intentionally conservative and named
`AutoRefine/default`:

- five macrocycles;
- reciprocal- and real-space coordinate refinement;
- group ADPs and eligible occupancies;
- automatic target selection, optimized XYZ/geometry and ADP weights, and the
  `n_gaussian` scattering table;
- no rigid body, individual ADP, TLS, NCS, reference-model, water-update,
  hydrogen-addition, or simulated-annealing strategy; and
- every logical processor available to the process.

When PostMR has no anomalous candidate, AutoRefine selects supported mean
observations. With a candidate, it requires an amplitude quartet
`F(+),SIGF(+),F(-),SIGF(-)`, forces anomalous observations to remain separate,
adds `group_anomalous`, and defines one exact atom selection per known model
site. A validated AutoSol MTZ supplies its discovered HL quartet at the same
time while the target remains explicitly `auto`. If AutoSol is approved but
F+/F- amplitudes are missing, a non-anomalous mean-data fallback may finish,
but the command returns a final error and marks the result for review.

Subprocess output is redirected to `phenix.refine.log`; the console receives
only a start notification and a curated result card. Rwork/Rfree history,
clashscore, bond/angle RMSDs, labels, selections, commands, outputs, and
acceptance decisions are written to the round report and a small metrics TSV.
Before computation, the same command is run through Phenix's dry-run parser in
the isolated round directory. An invalid parameter or input therefore produces
a small failed checkpoint and never launches the expensive refinement.

## Checkpoint graph

`AutoRefine/checkpoints.json` is the run-local registry. Nodes are immutable
and contain a stable ID, parent, recipe, status, model checksum, frozen input
references, metrics, compatibility result, and outputs. `SUCCESS` nodes become
current automatically. `REVIEW` nodes remain reusable but require explicit
selection; `FAILED` nodes remain visible and cannot be selected.

Checkpoint listing is read-only: it validates an existing
registry without rewriting it, or derives a validated PostMR root in memory
without creating files. Modifying operations still initialize the registry.

### Future graphical checkpoint-tree view

The immutable checkpoint registry is also the natural data model for a future
generic NASolve GUI. The GUI must visualize the existing graph rather than
create a second lineage system.

A **tree/graph mode** should show stage and model lineage across MR, PostMR,
AutoSol, ordinary refinements, Refine Doctor branches, manual imports, topology
surgery children and later re-MR attempts started from refined models. A
forty-refinement run should remain navigable as a graph rather than a flat
chronological list.

Useful graph operations are read/select/branch operations already compatible
with checkpoint semantics:

- pan/zoom through the lineage tree;
- collapse repetitive refinement chains into expandable groups;
- pin/bookmark important models as visible panels/cards;
- compare siblings or parent/child metrics without changing the current pointer;
- select any eligible node as the current checkpoint through the existing
  audited selection operation;
- launch `show`/Coot for one node without selecting it;
- branch a new refinement, Doctor, topology-surgery or later MR attempt from an
  explicit source node; and
- traverse back to an earlier model without deleting descendants.

The same graph component should be usable standalone in normal NASolve mode and
embedded in specialized workspaces such as Topo Net. Specialized workspaces may
add domain-specific overlays, but all state-changing actions still create
ordinary immutable checkpoint/model-lineage nodes.

The durable application-shell, scientific hierarchy, color/status semantics,
CLI-to-GUI coverage and capability-driven extensibility rules are defined in
[the GUI design contract](gui.md).

`nasolve checkpoints add RUN --name NAME` creates a bookmark without copying
the current node. Supplying `--model` imports a manual model into the run as a
new review child. `--mtz` is a deliberate observation replacement for an
external refinement and remains subject to label and Free-R validation before
the next round. `nasolve checkpoints use` changes only the current pointer;
all branches and provenance remain intact.

Automatic readiness requires a successful Phenix exit, a compatible model,
`Rwork < Rfree`, and `Rwork < 0.30`. These numerical tests permit continued
automation but do not replace visual inspection or user approval.

## Refine Doctor contract

Refine Doctor is a finisher and diagnostic layer, not an unconstrained
R-factor optimizer. It starts from the current or explicitly named checkpoint,
records that pointer, and requires every trial to be an immutable sibling of
the same source. A successful trial is recommended for inspection but is never
selected automatically by the engine. In an interactive CLI session the user
receives a final `[y/N/i]` prompt, with each answer confirmed by Enter. `i`
opens the exact recommended checkpoint through the ordinary Coot view path,
then repeats the prompt without waiting for Coot to close or rerunning Doctor.
Only `y` performs the ordinary audited checkpoint selection. `n`, blank input,
EOF or Ctrl-C at the prompt finishes with the pointer unchanged and prints
commands for inspection, later selection, and returning to the diagnosed
source. A Coot discovery or view error preserves the prompt for another choice.
Noninteractive sessions print commands without prompting or launching Coot.

Before launching Phenix it audits the authoritative Free-R array with the
Phenix/CCTBX runtime. The audit measures the independent Friedel-group count,
test fraction, mate consistency, and test-set population across resolution.
Inconsistent mate flags, fewer than ten independent free groups, or a guarded
fraction outside 2-20% stop all trials. A small or shell-sparse test set is
classified as noisy rather than invalid. Refine Doctor never tries alternate
test-flag values and never regenerates flags merely to reverse an R-factor
ordering.

An unavailable audit also stops all trials. A source that already passed the
numerical gate needs no default trials. Otherwise the declared sequence stops
at its first numerical pass, technical failure, or trial budget. Both
`--max-trials` (default five) and macrocycles per trial (default three) accept
1-10. Every attempted/skipped recipe and its reason is recorded.

Source anomalous benchmarks retain their wavelength, occupancy, B factor and
checkpoint. For selected anomalous observations, default Doctor trials use
explicit Henke f'/f'' values calculated with Phenix/CCTBX at the wavelength
in the MTZ header of the refinement observations (the dataset that holds
reflection data, read from the same `phenix.mtz.dump` output that selects the
labels). A recorded AutoSol wavelength is the fallback when the header has
none; a disagreement above 1e-4 A is recorded and the MTZ value is used. The
same wavelength gives the calculated f'' in every AutoRefine round's anomalous
diagnostics. The sequence tries ML and eligible MLHL with those values
fixed, then f''-only refinement, coordinate-only, and B-only fixed-scattering
fallbacks. Missing wavelength or failed calculation skips these recipes for
inspection; zero is not substituted. A successful f''-only execution can supply
a missing benchmark. Fixed values are distinguished from fitted values in the
metrics. The ordinary AutoRefine default remains unchanged.

The ordinary mean-data branches are:

- ML without HL phases, reciprocal-space XYZ, residue-group ADPs, no occupancy
  search, and anomalous-parameter refinement off; and
- the same parameterization with approved HL phases, allowing Phenix's
  automatic target to use MLHL.

An individual-ADP branch is eligible only at 3.2 A or better and with at least
three independent observations per modeled atom. This prevents low-resolution
DNA models from gaining a nominal R-factor improvement through an unsupported
parameter count. Trial definitions are declarative `RefineDoctorTrial`
records, so future project presets may supply additional reviewed strategies
without changing checkpoint or audit behavior. Coordinate-only and group-B-only
fallbacks follow; these are independent sibling trials, not chained stages.
Weight optimization is enabled only for the coordinate/ADP strategies being
refined, and occupancies stay fixed in every default Doctor recipe.

A branch is strictly successful only after a successful compatible AutoRefine
execution with finite valid R values, `Rwork < Rfree` and `Rwork < 0.30`.
Failed or incompatible outputs cannot win because they printed favorable
statistics. If none passes, Doctor returns `REFINE_DOCTOR_REVIEW` (exit 2),
leaves `recommended_checkpoint` empty, and separately identifies a usable
`inspection_checkpoint` by Rfree then Rwork. This descriptive ranking does not
quantify uncertainty or establish superiority. A small inversion no longer
endorses the original source as good enough. The additive `triage` and
`ranking` report fields explain budgets, skipped recipes, the stop and next
inspection steps. Campaign integration and diagnostic-driven restraint/data alternatives
remain planned in `campaigns.md`; immediate implementation state is tracked in
`development-handoff.md`.

## Stage-aware Coot views

`nasolve show` resolves the active workspace; an explicit run remains supported.
When a checkpoint registry exists, the selected current model takes precedence
over newer unselected attempts, including when current is manual or PostMR.
`nasolve show last DATASET` first chooses the highest numbered run. Explicit
`--stage` or `--checkpoint` opens a comparison without changing the current
pointer. Without a registry, the default prefers accepted AutoSol with a usable
density map, then completed PostMR, then Phaser. Profiles differ:

- AutoMR: Phaser model and Phaser MTZ;
- Metal (`--stage metal` on an AutoMR run, or a `DATASET/AutoSolMetal/run_NNN`
  directory opened directly): Phaser model, that search's (or the latest found
  search's) `overall_best_ha_pdb.pdb`, and the Phaser MTZ;
- PostMR: ReadySet model and dictionaries with accepted AutoSol density when
  available, otherwise an explicitly labelled Phaser map;
- AutoSol: ReadySet model and dictionaries, separate AutoSol HA model, and
  density-modified map coefficients; and
- AutoRefine: current refined checkpoint model, map coefficients, and its
  frozen dictionaries.

Manual checkpoints keep their own model and dictionaries. Maps may be inherited
through the checkpoint ancestry only when frozen observations agree; the source
label states they were not recalculated for the manual coordinates. Raw HL-only
phase data are never substituted for a density map. Declared sources and their
checksums take precedence; ambiguity, corruption, or a missing declared map is
an error. Legacy AutoSol runs retain bounded discovery within their own stage.

Every view retains ordinary `--auto` map loading and adds a run-local Python
startup script. In Coot, `valid_labels` checks the exact `ANOM`/`PHANOM` pair.
The script uses `make_and_draw_map` with no weights and the difference-map flag,
then names the overlay and sets its initial contour to 3 sigma. It adds no
phase shift. Existing maps with the same MTZ and column pair are reused, and
the original refinement-map and scroll-wheel selections are restored.

Candidates are the selected map MTZ followed, for refinement checkpoints, by
the explicitly recorded `refinement_reflections` MTZ if different. Both follow
the existing relocation and checksum rules; a corrupt or missing declared
source stops launch. Manual checkpoints inherit these candidates from the
same observation-compatible ancestor as their ordinary maps. Other stages
check only their selected map. There is no search across refinement rounds,
dependence on an anomalous-refinement flag, or substitution of HL coefficients.
Missing columns simply leave the normal view in place. Optional Coot API/read/
display failures are logged without closing the ordinary view.

`launch.json` has an additive `anomalous_map` request with paths, labels,
contour, script path, and `status: check-in-coot`; this is not a claim that an
asynchronously started GUI loaded the map. `NASOLVE_ANOMALOUS_MAP` JSON lines
in `coot_gui.log` record `loaded`, `reused`, `absent`, or `error`. No checkpoint
or scientific report schema changes, regeneration, or dependency installation
are required. Tests execute the generated script against a simulated Coot
Python API; graphical validation on the user's Coot 1.3.3 / Phenix 2.2.1-6174
outputs remains a live check, not a claimed execution in the test environment.

Every graphical process starts in `RUN/CootGUI/STAGE/` with its backup
directory redirected underneath the same pen, plus a checkpoint subdirectory
when selected. Console output and `launch.json` identify the run, checkpoint,
model, and map source. Run discovery uses numbered directories and reports,
never filesystem modification times.

## Campaign planning contract

`presets.load_preset` validates a schema-1 TOML policy and its declared resource
bytes. The packaged `5w6w` preset covers the existing standard workflow. Unknown
settings and unsupported policies stop loading. Paths resolve within the preset
root, and exact source and resource hashes accompany the canonical configuration
fingerprint. Python 3.10 uses the conditional `tomli` dependency; newer Python
uses `tomllib`.

The optional `[chemistry] terminal_phosphate_sites` list is a typed, explicit
site declaration. The W frame points to the packaged `5w6w` card (now version
1.1.0) for standalone input resolution. Campaign planning passes its selected
card to the same resolver; a custom card does not silently inherit the built-in
card's chemistry. The card is the single source of the D:1 declaration.

New resolutions freeze `phosphate_intent` (schema 1, effective list, source
`recipe`/`dataset`/`none`, recipe id/version/frame and source/config SHA-256)
in each campaign `effective_config` and AutoMR `post_mr_plan`. The existing
`allow_op3_sites` field is the operational list, so PostMR, ReadySet validation,
AutoRefine, Doctor and manual children use the same permission unchanged.
Campaign loading validates recipe provenance against the frozen preset and
checks the effective list against its source. Execution does not consult the
current card to add permissions. CLI plan/status reports sites and their source.
`format_intent` always emits even an empty `allow_op3_sites` field, preventing
re-read snapshots from acquiring later defaults.

Old runs/plans lacking this additive provenance still use their original frozen
list (or no consent if absent). Do not patch an existing run or campaign in
place just because a recipe changed; prepare a fresh attempt with explicitly
selected chemistry. Failed ED run_002 remains historical evidence.

`campaigns.plan_campaign` composes existing AutoMR input discovery, intent parsing,
alias resolution and catalogue selection without allocating a scientific run.
Discovery is limited to immediate dataset children. A malformed or incomplete
dataset becomes `BLOCKED` while other entries are retained. The planner freezes
the dataset list, explicitly resolved settings, checksums and duplicate MTZ
groups. Per-file references use a campaign anchor and relative path; external
selected models and preset resources are snapshotted into plan storage.

`NASolveCampaign/plan.json` is an immutable schema-1 planning record. Publication
must not overwrite an existing plan or incomplete planning directory. The
manifest is exposed only after its resource snapshots are complete. Status
validates the existing record and input identities without rewriting it or
admitting new datasets. Relocation tests remove the old campaign, source preset
and model catalogue before checking the copied plan.

Planning commands do not discover Phenix/Coot, touch workspace settings, select
checkpoints or execute stages. `DISCOVERED` confirms only input and model
selection. Full model/symmetry/array validation and downstream chemistry gates
remain with the existing engines. The executor verifies the frozen plan before
consuming it and stores execution progress separately. See
[campaign planning](campaign-planning.md).

## Campaign execution contract

`campaign_execution` composes one standard W/5W6W path per dataset through the
existing scientific engines. The foreground executor is local and sequential;
individual Phenix stages retain their normal processor allocation. It consumes
the plan's resolved configuration and resource snapshots rather than choosing
a new model from a possibly changed source catalogue.

The immutable plan remains schema 1. Separate schema-1 execution state under
`NASolveCampaign/execution/` records the plan fingerprint, exact attempt/run
ownership, stage progress, checkpoint, diagnostic and inspection requirement.
Atomic updates, campaign/dataset locks, process identity and heartbeat records
protect against competing execution and expose interrupted work. Stage workers
execute in owned process groups so cancellation can terminate their external
tool children. A pause request finishes the active stage before stopping.

Resumption verifies saved stage results and checksummed artifacts before
continuing. It never adopts an unrelated newest run, guesses artifacts by
basename, or overwrites an incomplete stage. An explicit retry preserves the
old attempt and schedules a new preflight with a fresh numbered run. Live or
uncertain prior processes prevent duplicate work. Whole-campaign relocation
requires all workers to be stopped first and preserves relative artifact
references; tool configuration is rediscovered on the current computer.

MR review, unaccepted AutoSol phases and refinement review stop the affected
dataset for inspection while other selected datasets continue. A numerical
pass records its selected checkpoint as `SOLVED` with both
`numerical_success = true` and `inspection_required = true`. Final structural
approval, automatic Doctor selection and deposition remain separate future
layers. The existing observation, Free-R, chemistry and checkpoint selection
gates are preserved. See [campaign execution](campaign-execution.md) for the
command and recovery contract.

## Project preset direction

The first planning preset makes a bounded set of project settings declarative.
Further frame and project policy should use data rather than additional
conditionals keyed to names such as `5W6W`. A later preset manifest beside each
frame catalogue can declare:

- model providers, exact-pair catalogues, fallbacks, and copy/symmetry policy;
- standard sites, chain sequences, and restraint resources or modes;
- reviewed construct-registration/cut recipes;
- metalation recipes and anomalous-element policy;
- the AutoSol sequence resource and phasing defaults; and
- imported external-provider provenance/capabilities.

The orchestration layers consume a frozen model plus declared capabilities.
This permits a new experimental campaign to ship a versioned preset directory
without changing common run allocation, provenance, safety gates, Coot/Phenix
isolation, or downstream reporting.

Model generation is not part of NASolve's runtime contract. A separate optional
NAPrep package may manage design/data records and externally generated model
candidates (including AlphaFold outputs), but NASolve consumes only the frozen
handoff and then owns the downstream crystallographic/campaign decision tree.

## Next development priorities

Earlier end-to-end local validation is preserved in
`history/validation-dohu.md`. Current validation state and the immediate
implementation edge are maintained in `development-handoff.md`.

The campaign planner freezes inputs and project policy; the sequential executor
adds guarded stage composition and saved progress. The read-only Campaign
Doctor provenance/comparison prerequisites are now implemented.

The next structural robustness layer is construct registration: a cheap
non-mutating Registration Scout before Phaser, authoritative ASU Registration
on the MR solution before PostMR, a guided Registration Net for ambiguous
cases, and bounded Registration/Recut Rescue candidates only when ordinary MR
or representation mapping needs them. Current clean W runs must retain an
identity-like fast path.

Reviewed Campaign Doctor eligibility/rescue policy and richer inspection
summaries remain the next major orchestration steps after this registration
foundation.

Campaign orchestration should reuse frozen inputs, immutable numbered runs,
and checkpoint lineage, with resumable per-dataset progress and explicit
reporting of a blocked dataset. Supported residue mappings should run without
per-residue prompts; unavailable or ambiguous chemistry must retain a useful
diagnostic. Numerical readiness, recommended branches, and user inspection
must remain distinct states. Project presets should carry scientific policy
instead of adding project-name branches to the orchestration code.

NARestraints follow-up is separate from the completed PostMR connectivity fix.
User inspection of Q:E found a collapsed modified-base wobble arrangement with
suspect angles despite plausible contact distances; sulfur–nitrogen contact
targets also need review. Resolve these against atom identities, maps, and the
generated restraints before changing angle or distance targets. No such
restraint-target change is included here.

Reports should retain raw validation scores together with model and lattice
context. Simon identifies the DOHU structure as an intentionally interconnected
symmetry-related crystal, where the reported clashes are expected. That
interpretation belongs alongside inspection evidence; it does not require
altering the recorded clashscore. Future campaign acceptance policies should
be preset-specific and should preserve the distinction between numerical
acceptance and structural assessment.
