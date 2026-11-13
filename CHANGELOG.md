# Changelog

All notable user-visible, compatibility, schema, and reproducibility changes
to NASolve are recorded here.

The format follows Keep a Changelog. Until the first tagged release, historical
entries are reconstructed from repository history.

## [Unreleased]

### Fixed

- The metal search (`nasolve autosol RUN --metal ELEMENT`) gave AutoSol the
  plain Phaser model, so heavy atoms that belong to the requested bases were
  found as metal sites: in iCiU_Ag the iodines of 5-iodo-dC and 5-iodo-dU
  appeared as two extra "Ag" sites beside the real one. The residues the
  AutoMR run asks for (standard pair, explicit mutations or sequences) are now
  built into a copy of the Phaser model first, with PostMR's own mutation step
  (headless Coot, coordinates restored), in
  `AutoSolMetal/run_NNN/PartialModel/partial_model.pdb`, and that copy is
  AutoSol's partial model. Coot is required only when something must be built;
  `--original-model` keeps the old behaviour. PostMR's mutation step moved into
  `postmr.plan_model_mutations` / `apply_model_mutations` without changing
  PostMR. Schema: the metal-search report gains `inputs.partial_model` and
  `partial_model` (source, mutations, Coot log, coordinate restoration); site
  distances, the Coot view and the merge base use the partial model, and the
  merged model keeps its Phaser-based name. Coot and AutoSol simulated in
  tests; not yet run on real data.
- The metal search (`nasolve autosol RUN --metal ELEMENT`) stopped with
  "AutoSol needs MR_frames/FRAME/seq_base.txt or complete run sequences" for
  any AutoMR run without a frozen frame sequence, such as every 3GBI-frame run
  (only the W catalogue ships `seq_base.txt`), because PostMR sequences do not
  exist yet at that stage. It now falls back to the nucleotide chains of the
  Phaser model, which hold every strand in the asymmetric unit (residues with
  C1'; curated modified residues as their parent; DU as T; others left out with
  a warning), and `--sequence FILE` gives a sequence explicitly (metal search
  only). A frozen sequence that is missing or fails its checksum is still an
  error. Schema: metal-search `report.json` gains `sequence_warnings`;
  `sequence_source` can now be `--sequence NAME` or `derived from the Phaser
  model (...)`. The command prints the source. Phenix simulated in tests.
- AutoMR's `phaser.eff` quoted paths with JSON escaping, so non-ASCII
  characters were written as `\uXXXX`. PHIL keeps such sequences literally,
  and Phaser stopped with "Couldn't find the file" for any dataset whose path
  contained, for example, Chinese characters. Paths are now written as UTF-8
  and only `\` and `"` are escaped, matching `libtbx.phil` quoting (as
  AutoRefine already did). Checked against the Phenix 2.2.1-6174 PHIL parser
  (parse only, not a Phaser run). Existing failed runs are unchanged; rerun
  AutoMR to create a new run.
- AutoRefine and Refine Doctor took the X-ray wavelength only from an AutoSol
  record, so a structure refined without AutoSol (for example 7CU_Ag, whose Ag
  site came from the metal search) reported no calculated f'' and Doctor
  skipped every fixed and f''-only scattering recipe. The wavelength now comes
  from the MTZ header of the refinement observations, read from the same
  `phenix.mtz.dump` output that selects the labels (the dataset that holds
  reflection data, as for AutoSol). A recorded AutoSol wavelength is used only
  when the header has no usable one, and a disagreement of more than 1e-4 A is
  recorded as a warning, with the MTZ value used. Schema: each AutoRefine
  round's `inputs` gains `wavelength` and `wavelength_source` (`method`,
  `mtz_dataset`, `mtz_error`, `autosol_wavelength`, `warnings`); anomalous
  diagnostics and calculated scattering groups record `wavelength_source` as
  the method used. Existing rounds are unchanged. Phenix programs simulated in
  tests; the parser is the one that read 1.00743 A from the 7CU_Ag MTZ in its
  metal search (Phenix 2.2.1-6174).

### Added

- `nasolve autofdp MODEL.pdb DATA.mtz` (steps 1-3 of the f'' occupancy
  calibration specified in `docs/autofdp.md`, after the mmDNA supplements:
  Adv. Mater. 2023, JACS 2023, Small 2025, Sci. Adv. 2026). It needs only a
  starting model and an MTZ; it reads no run report, AutoRefine round or
  checkpoint. From the MTZ (`phenix.mtz.dump`): the anomalous observations
  (F(+)/SIGF(+)/F(-)/SIGF(-) or the intensity equivalent; required), the
  Free-R label (never generated; test value 1 for `R-free-flags`, 0 for the
  CCP4 labels) and the wavelength (MTZ header, or `--wavelength`). From the
  model: the metal HETATM sites, found as in PostMR.
  - Step 1 (B start): the model is copied unchanged to `source_model.pdb` and,
    with each metal's B set to the atom-weighted mean B of its neighbouring
    nucleotides (NARestraints metal-pair search: the bridged pair; otherwise
    base contacts within 3.6 A; otherwise nucleotides within 6 A), to
    `start_model.pdb`. Only the B columns change. `--site`, `--b-start
    model|VALUE`.
  - Step 2 (f'' round): first `phenix.ready_set` on `start_model.pdb` with
    PostMR's actions (no hydrogens, dictionaries for unknown residues, no
    geometry optimisation), given NASolve's curated dictionaries for curated
    residues in the model and any `--cif` files; its model must keep every
    atom, coordinate, occupancy and B and becomes `readyset_model.pdb`, and its
    dictionaries (PostMR's precedence, `effective_restraints`) go to
    phenix.refine. Then one phenix.refine round on `readyset_model.pdb` with
    every strategy off except `group_anomalous` (no coordinate, real-space, B
    or occupancy refinement; one group per metal site, f' fixed at its cctbx
    Henke value for the wavelength, f'' refined from 0; `--cycles`, default
    3), with parameters from AutoRefine's writers and a `--dry_run` check
    first. AutoFDP prints each refined f'' against theory (deviation, within
    `--tolerance`, default 1 e-), checks that coordinates, occupancy and B did
    not move, and opens the round's model and maps in Coot (`--no-coot`,
    `--coot`). `--theory ELEMENT=FDP` replaces an element's theoretical f''.
    `--cif` adds dictionaries for ReadySet and phenix.refine. Schema:
    `report.json` gains `readyset` and `steps.readyset`.
  - Step 3 (occupancy adjustment): the calibrated sites one at a time, the one
    furthest from theory first. Occupancy goes down when the refined f'' is
    below theory and up when above (`propose_occupancy`: proportional first
    step q x f''/theory, then a secant on 1/q through the site's last two
    occupancies; at most 4-fold and at least 0.01 per step, within
    `--occupancy-floor`, default 0.05, and 1.0), written into a copy of the
    previous round's model (`round_NNN/input_model.pdb`, only that occupancy
    changed), and the same f''-only round runs again (same parameters,
    cycles and dictionaries; f'' from 0). A site stops when it is within
    `--tolerance`, after `--max-iterations` adjustments (now default 5, per
    site; 0 runs round_001 only), or at full occupancy / the floor with f''
    still on the wrong side. The last round's model is `final_model.pdb` and
    opens in Coot with its maps; the command prints every adjustment and a
    final table. Status `AUTOFDP_CONVERGED` (exit 0) when every calibrated
    site is within tolerance, `AUTOFDP_REVIEW` (exit 2) with per-site reasons
    otherwise, `AUTOFDP_FAILED` (exit 1) when a round fails. Schema (additive,
    still 2): each `rounds` entry gains `input_model`, `adjusted_site`,
    `occupancy_before`/`after`, `method` and `note`, and each site its
    `alternate`; `report.json` gains `occupancy` (method, limits, every
    adjustment, per-site outcome and history), `final_model`, and
    `settings.max_iterations`/`occupancy_floor`; `steps.occupancy_adjustment`
    is now `done`, `review`, `failed` or `not run`. The `decide_site_update`
    stub is replaced by `propose_occupancy`.
  - Output: a new `run_NNN` in `AutoFDP/` beside the MTZ (`--output` for
    another folder). Schema: `report.json` (`workflow: "autofdp"`, schema 2)
    with portable, checksummed references to the inputs (dataset anchor) and
    outputs (run anchor), what was read from the MTZ, the theory, the B start,
    `rounds` and `occupancy`; each round's status is `AUTOFDP_FDP_ROUND_DONE`
    or `AUTOFDP_FAILED`.
    `AutoFDP` joins `AutoMR` and `AutoSolMetal` as a dataset-level run series
    for path anchoring. `--plan` does the preflight (running only
    `phenix.mtz.dump` and `phenix.python`) and writes nothing.
  - Not implemented yet: the identity check, B confirmation, the final round
    and the map-based choice between B and occupancy (documented stubs in
    `nasolve.autofdp`); `--element` is refused. phenix.refine, phenix.mtz.dump, phenix.python and Coot are
    simulated in tests; not yet run on real data.
- `nasolve autosol RUN --metal ELEMENT` searches for metal ions by MR-SAD
  AutoSol directly after Phaser, with the original Phaser model as partial
  model (`autosol.phasing.partpdb_rms`, default 1.0 A, `--rmsd`), extra metals
  via `autosol.phasing.mad_ha_add_list`, and no model building. Each search
  writes a new `DATASET/AutoSolMetal/run_NNN`, beside `AutoMR` and the input
  data (earlier development builds wrote `AutoMR/run_NNN/AutoSolMetal/
  attempt_NNN`; those schema-1 records stay viewable). Outcomes
  `METAL_SITES_FOUND`, `METAL_NOT_FOUND` (normal) or `METAL_AUTOSOL_FAILED`.
  Found sites are listed with their nearest residues and opened in Coot with
  the Phaser model and map (`--no-coot` to skip; reopen with `nasolve show
  DATASET/AutoSolMetal/run_NNN` or `nasolve show RUN --stage metal`). Schema:
  the search's own `report.json` (`workflow = "metal_autosol"`,
  schema_version 2, `source_run`, references anchored to the search directory
  or the dataset) and an additive `report.json["metal_autosol"]["searches"]`
  index in the AutoMR run; the run's stage/status are unchanged and sites are
  never added to the model. `run_context.dataset_directory` now also treats
  `AutoSolMetal/run_NNN` as a dataset-level run.
  `nasolve autosol` without `--metal` is unchanged. Tested with a simulated
  AutoSol only; live Phenix 2.2.1-6174 validation pending.
- PostMR writes metal ions as anomalous scatterers. After ReadySet it finds
  HETATM metals with atomic number >= 21 outside nucleotides (Ag, Hg, Zn,
  ...; not Na, Mg, K, Ca) using NARestraints `restraints.anomalous`
  (anomalous-scatterers branch), records them in
  `postmr.anomalous.metal_sites` (same fields as `candidates`, plus
  `kind = "metal"`), and writes `Restraints/anomalous_scatterers.phil` with
  NARestraints' `anomalous_phil_block` (one group per element, f'/f'' from 0
  and refined, `group_anomalous` in the strategy) for hand-run
  `phenix.refine`. AutoRefine now treats `metal_sites` like nucleotide
  heavy-atom candidates: it selects F(+)/F(-) data and refines one
  anomalous group per site. The PHIL is deliberately not in the AutoRefine
  restraint stack, to avoid duplicate groups. `autosol_required` still
  reflects nucleotide heavy atoms only. With an older NARestraints, PostMR
  stops with a clear message when such a metal is present. For PostMR runs
  made before this change (no `metal_sites` key), AutoRefine and the refine
  doctor scan the checksum-verified PostMR prepared model the same way, in
  memory, without rewriting the frozen PostMR record. Checked on the real
  7CU_Ag PostMR (Ag at B:101 after Coot's merge: `chain B and resid 101 and
  name AG`); phenix.refine itself not yet run.
- dU (`DU`, the DNA `U` token, for example `--pair C:U`) is now a curated
  PostMR residue. `data/ligands/DU.cif` is the official wwPDB CCD definition
  copied unchanged from Phenix 2.2.1-6174
  (`chem_data/chemical_components/d/data_DU.cif`, DNA LINKING, wwPDB
  modified date 2011-06-04). PostMR builds DU on its parent DT: Coot overlays
  DU on the DT already in the model and every DU atom takes the DT
  coordinates, so the result is the original dT without its methyl C7.
  ReadySet now receives DU.cif instead of generating a non-polymer DU
  dictionary. Previously PostMR stopped with "No local CCD dictionary for DU".
- PostMR now plans its mutations and checks every needed ligand dictionary
  (targets and curated residues already in the model) and the Coot
  executable before creating `PostMR/`. Previously a missing dictionary or
  Coot left a half-made `PostMR/` folder that made every later attempt fail
  with "refusing to overwrite".
- The metal Coot view (after `nasolve autosol --metal` finds sites, or
  `nasolve show DATASET/AutoSolMetal/run_NNN`) adds a **NASolve** toolbar
  button, "Merge metal sites and save", and the Coot Python command
  `nasolve_merge_metal_sites()`. It merges the open (possibly edited)
  `overall_best_ha_pdb.pdb` molecule into a copy of the starting Phaser model
  and writes `<model>-merged.pdb` (for example `mr_solution-merged.pdb`) in the
  source AutoMR run's `Phaser/` folder, adding `-2`, `-3`, ... rather than
  overwriting. Occupancies above 1 are set to 1.00 and listed; sites whose
  recorded symmetry-aware nearest residue is clearly closer than the site as
  placed are flagged as belonging next to a symmetry copy. The open molecules
  are not changed. The merge is recorded in the AutoMR run's
  `report.json["metal_merge"]` (schema_version 1: `models` with run-anchored,
  checksummed model, base Phaser model and metal-search references, sites
  merged, occupancy changes and warnings; `selected` = newest), written
  atomically from Coot. `launch.json` records `metal_merge`, and the Coot log
  records `NASOLVE_MERGED_MODEL`. Coot calls follow the Coot 1.3.3 Python API
  bundled with Phenix 2.2.1-6174; tested against a simulated Coot (not yet in
  a live Coot session).
- PostMR starts from the selected merged metal model when one is recorded
  (`nasolve postmr RUN --original-model` keeps the plain Phaser model). A
  selected model that is missing or changed since merging stops PostMR before
  `PostMR/` is created. The PostMR report adds `input_model_source`
  (`phaser` or `metal-merged`) and the `metal_merge` entry used; `input_model`
  and `input_sha256` name the merged file. With the ion in the model,
  NARestraints (anomalous-scatterers branch) finds the metal-mediated pair and
  skips its base-pair and stacking restraints; checked on the 7CU_Ag Phaser
  model plus its AutoSol Ag site (B:4 / A:12 via Ag). Runs that already have
  PostMR are unchanged. The metal searches still use the original Phaser
  model, so a merged ion never validates itself.
- `nasolve show PATH` looks a relative `PATH` up from the current folder and
  then from the active workspace dataset, so `nasolve show
  AutoSolMetal/run_001` works outside the dataset folder. When neither has a
  `report.json`, the error names both places and says which folders are
  expected, instead of the misleading "AutoMR report is missing".
- AutoMR now always tells Phaser which data to use: `-F` (merged
  amplitudes, the default) or `-I` (merged intensities, as the Phenix GUI
  uses); `--phaser-data F|I` is the long form. Or set `phaser_data = F|I` in
  `[automr]` of
  `nasolve.txt`; the command line overrides the file. Previously `phaser.eff`
  had no `labin` and `phenix.phaser` chose for itself (`F,SIGF` whenever the
  MTZ had amplitudes). The preflight lists the MTZ column types with
  `phenix.mtz.dump` and groups them as `iotbx.mtz` does: `I` takes the first
  merged intensity array (`J` then `Q`; `IMEAN,SIGIMEAN` for STARANISO files);
  `F` takes `F,SIGF` or `FOBS,SIGFOBS`, else the first merged non-`FC`
  amplitude array (phenix.phaser's own preference). The preflight stops before
  allocating a run if the requested columns are missing, and now needs
  `phenix.mtz.dump` in nonstandard mode too. The choice is always written to
  `nasolve.input.txt`. Schema: additive `report.json["phaser_data"]`
  (schema_version 2: `requested` F|I, `labin`) and
  `execution.phaser.reflection_data`, which also records Phaser's own `LABIN`
  echo; a run whose log reports the other kind of data is `MR_FAILED`. Reports
  frozen before this record, and schema-1 records from the first version of
  this option (`auto`/`intensities`), still run as they did. Campaigns use `F`
  and block a dataset that sets `phaser_data = I` rather than ignoring it.
  Column selection was checked against the column header of a real
  STARANISO MTZ and the Phenix 2.2.1-6174 `iotbx.mtz` and `phenix.phaser`
  sources; Phaser itself was simulated in tests (no live run yet).
- `nasolve automr ... --execute --show` opens Coot when Phaser has written a
  model and map, whatever the TFZ result: Phaser's model and its map
  (`Phaser/mr_solution.pdb` and `.mtz`, byte-identical copies of
  `PHASER.1.pdb`/`.mtz`) in an isolated `CootGUI/automr` session, as
  `nasolve show RUN --stage automr` does; `--coot` overrides the executable.
  Without `--show` nothing opens (unchanged default) and AutoMR prints the
  `nasolve show` command instead. A Coot problem prints an error but never
  changes the AutoMR status or exit code. Tested with a simulated Coot launch
  only.

### Changed

- The MTZ header is now the source of the space group and the X-ray
  wavelength, and `Data_1*.cif` and `summary.html` are optional. Datasets
  without them (for example the 2021 iCiU_Ag autoPROC run, which has only the
  MTZ files and `report.pdf`) are accepted. When either file is present it is
  still a cross-check and must agree with the MTZ. Standard-frame symmetry
  compares the `phenix.mtz.dump` symbol and symmetry matrices with any CIF and
  summary entries; AutoSol and the metal search read the wavelength of the MTZ
  dataset that holds reflection data (the `HKL_base` placeholder is ignored),
  and a present `summary.html` must agree within 1e-4 A (one without a
  wavelength is recorded as not cross-checked). The MTZ stores the wavelength
  to 5 decimals, so `autosol.lambda` may differ from the summary's 6-decimal
  value by up to 5e-6 A. Schema: `report.json["inputs"]["metadata"]` and
  `["summary"]` may be null; `symmetry.evidence` gains `compared_sources`
  (`cif_symbols` empty, `summary_symbol` null when absent); the AutoSol and
  metal-search `wavelength_source` is now an object (`method`, `path` of the
  MTZ, `mtz_dataset`, `summary_cross_check`). Campaign plans omit absent
  optional inputs. Existing runs are unchanged. Supersedes the unmerged
  `wip/cif-wavelength` branch (wavelength from the CIF). Checked against the
  Phenix 2.2.1-6174 `iotbx.mtz` dump format and the iCiU_Ag and FT_Hg MTZ
  headers; Phenix programs simulated in tests.

- Added schema-1 semantic `terminal_geometry_protection` checkpoint provenance. Protected sites, sigma, Phenix-derived ideal source, mechanism, and the run-local protection restraint artifact now inherit with refinement/manual children; attempts to stack a second protection record fail closed.
- Generalized the OP3-specific policy into an explicit standard-phosphodiester backbone contract. `five_prime_phosphate_sites` is the preferred user-facing name (legacy `allow_op3_sites` remains readable). PostMR now treats a requested 5'-terminal phosphate as the complete P/OP1/OP2/OP3 group, preserving a complete group, completing missing OP3 from existing P/OP1/OP2, or seeding a whole missing group from O5'-C5' with recorded idealized starting geometry. Partial ambiguous groups still fail closed.
- Terminal-phosphate native geometry remains owned by Phenix. NASolve now audits declared terminal phosphates from Phenix geometry output and provides provenance-tracked local `action = change` protection using Phenix-derived ideals; NARestraints remains responsible for reviewed pairing/stacking geometry and does not duplicate native terminal-phosphate restraints. Automatic first-AutoRefine protection is now implemented and live-validated: NASolve obtains native ideals from a pre-refinement `phenix.pdb_interpretation` geometry snapshot, applies six local sigma=1-degree P-centered `action = change` restraints, freezes the source geometry/log/protection artifact into checkpoint provenance, and inherits that protection without duplicate stacking.
- Added site-scoped `experimental_passthrough` for explicitly declared non-standard backbones. Standard phosphate rules are skipped only at those sites; no custom linkage is inferred. Passthrough requires explicit user authorization, persists in run/campaign provenance, and has a two-step `backbone-review` Coot/confirmation workflow that records the inspected model hash without erasing provenance. Reviewed custom backbone recipes and 3'-phosphate construction remain future work.

- Added explicit `[chemistry] terminal_phosphate_sites` to recipe cards. The
  built-in W/5W6W card is now version 1.1.0 and declares D:1, confirmed by Simon
  as a designed 5-prime phosphate. Standalone W selection and campaign planning
  resolve this same card; custom cards declare their own lists. Dataset
  `allow_op3_sites` replaces recipe sites, including an explicit empty override.
- Freeze recipe phosphate origin/id/version/hashes beside the effective sites in
  AutoMR and campaign records. Report sites in AutoMR and campaign plan/status.
  Existing frozen runs/plans never gain permissions from a later recipe version.
  Fresh PostMR attempts may construct an explicitly requested 5-prime terminal
  phosphate under the current guarded chemistry contract; old runs remain immutable.

- OP3/O3P terminal intent is strictly explicit and is frozen into
  AutoMR/campaign provenance. For a requested 5-prime terminal site, PostMR may
  preserve a complete group, add missing OP3 to P/OP1/OP2, or construct a wholly
  missing P/OP1/OP2/OP3 group from the local sugar frame. Ambiguous partial
  groups, unrequested unlinked OP3, cyclic/incoming-link conflicts, and an extra
  oxygen on an internal phosphate still fail closed.
- Replaced the bundled raw 1AP CCD graph with the tested parameterized monomer
  library dictionary, preserving numerical values and adapting its group to
  DNA. New PostMR 1AP profiles generate residue-selected OP3 modifications and
  freeze checksummed refinement/view artifacts. ReadySet cannot supersede the
  reviewed 1AP definition or swallow modification files during checkpoint
  creation. Normal dictionaries and coordinate/raw run evidence remain intact.
- Added phosphate/profile gates to AutoRefine, Doctor triage and manual imports;
  retained inherited modifiers after relocation, with missing or changed
  artifacts failing closed. Combined dictionary inputs merge their component
  lists instead of concatenating conflicting data blocks.
- NASolve now installs the released NARestraints `v1.1.2` dependency. No DE,
  sulfur pair-target, observations, Free-R, or NARestraints workbook changes are
  implied by this backbone work. ReadySet/Phenix interpretation and live
  refinement of a wholly constructed terminal phosphate have now been exercised;
  the earlier 1AP integration diary is archived under `docs/history/`.

### Added

- Added an explicit sequence-reference target compiler with the reviewed
  W-scaffold reference, chain/residue correspondence, thread/dataset overlays,
  assignment provenance and read-only inventory comparison.
- Added opt-in standalone `[automr] sequence_reference` integration. Fresh
  AutoMR runs freeze exact reference bytes and effective targets; PostMR reuses
  existing mutation routes and audits complete prepared residue identities.
  Snapshot/intent drift fails closed. Ordinary W defaults and legacy runs are
  unchanged. Campaign planning now snapshots an explicitly selected per-dataset
  reference into content-addressed resources; campaign preflight consumes that
  frozen copy rather than reopening the installed or dataset source.
- Added explicit campaign sequence threads through root
  `nasolve-campaign.toml`. Threads freeze named dataset membership, inherit the
  reviewed W sequence reference plus shared chain/site overlays, and preserve
  provenance into run-local sequence-family targets. Dataset declarations remain
  more specific. Membership does not imply MR-model compatibility or reuse.
- Added explicit standard-frame model providers through `[automr] model`
  and `automr --model`. A forced PDB may come from the dataset or selected frame
  catalogue, is frozen with structured provider provenance, and retains the
  frame's symmetry/recipe/sequence context. Forced filenames never imply pair
  compatibility; required W PostMR sites are checked literally before MR.
  Campaign plans snapshot the exact model bytes and provider so execution remains
  portable after the original model/catalogue disappears.
- Added checksum-verified search-model mismatch provenance for runs with a
  complete explicit residue target. The frozen comparison records chain/site
  inventory, missing/unexpected sites, literal residue mismatches, and the target
  source layer for each expected PostMR correction. It is descriptive only:
  no compatibility score, mutation-route guarantee, donor eligibility, or
  automatic model reuse is inferred.
- Added checksum-verified model compatibility fact sheets for every fresh
  AutoMR run. They record provider/frame declarations, exact candidate
  chain/residue inventory, target-dependent site/identity relations when
  available, mirror-transform state, recipient terminal/backbone chemistry and
  recipient symmetry/copy-number evidence. Unknown dimensions remain explicitly
  unknown. Schema 1 forbids scores, overall compatibility verdicts, donor
  eligibility and automatic reuse authorization.
- Added explicit provider-side `[automr] model_family` declarations for named
  search models. The declaration is bound to the selected model, frozen through
  campaign relocation, and compared literally with an explicit recipient
  sequence-reference identifier when both exist. Frames, filenames, sequence
  threads and similarity never create family membership implicitly, and the
  resulting SAME/DIFFERENT fact still carries no donor-eligibility decision.
- Added read-only checkpoint model-candidate descriptors as a prerequisite for
  future Campaign Doctor rescue. Descriptors verify checkpoint lineage, model
  and source-observation checksums, literal residue inventory, frozen chemistry,
  source-provider provenance and complete sequence-family target context when
  available. Checkpoint models are compared afresh with the frozen target.
  Local checkpoint reusability is recorded separately from cross-dataset donor
  eligibility, which remains null; no model is selected, exported or reused.
- Added a read-only donor-checkpoint versus recipient-run comparison layer.
  It re-verifies Fern donor provenance, consumes only the recipient run's frozen
  AutoMR target/context, compares literal donor residue inventory to the
  recipient target, and records mode/frame/reference/mirror/symmetry context as
  separate descriptive dimensions. Scores, ranking, donor eligibility,
  recipient compatibility, rescue authorization and automatic reuse remain
  explicitly unset/false.
- Clarified that experimental-family membership and MR-model reuse eligibility
  are separate, and that a metal-adapted scaffold does not imply metal pairing.

- `show` now opens an additional anomalous difference map at 3 sigma when
  Coot finds `ANOM`/`PHANOM` in the selected map MTZ or that checkpoint's
  declared refinement MTZ. This also applies to Refine Doctor inspection and
  compatible manual-checkpoint ancestors, irrespective of whether f'' was
  refined. Existing overlays are reused; ordinary maps and their controls are
  preserved. Startup diagnostics record loaded, reused, absent, or error.
- Refine Doctor now tries a bounded mean-data recovery list including
  coordinate-only and group-B-only siblings. Anomalous trials explicitly fix
  wavelength-calculated scattering or refine f'' alone. `--max-trials` defaults
  to five; it and `--cycles` accept 1-10. Trials stop at the first numerical
  pass, technical failure, or budget. Additive reports preserve scattering
  provenance, skipped-recipe reasons, descriptive ranking and inspection steps.
- Unresolved Doctor inversions now remain `REFINE_DOCTOR_REVIEW` (exit 2), with
  a separate inspection candidate. They no longer endorse the source or choose
  a review branch by its R-factor gap. Unavailable Free-R audits stop trials;
  failed/incompatible outputs cannot win by favorable printed statistics.
- Documented the proposed Refine Doctor triage sequence for noisy mean-data
  refinement, conditional anomalous recipes, and campaign stop/inspection
  rules, distinguishing the implemented standalone options from future work.
- Added foreground sequential `campaign run`, stage boundaries, pause requests
  and explicit fresh-attempt retries for frozen schema-1 plans. Scientific
  stages retain their existing gates; blocked and review datasets stop locally
  while other selected datasets continue. Separate schema-1 execution records
  preserve attempt/run ownership, stage results, process identity and recovery
  diagnostics. Resume verifies completed artifacts, refuses to overwrite
  incomplete work, and preserves checkpoint lineage. Campaign status now shows
  exact runs, next stages and checkpoints; numerical success still requires
  model/map inspection. Execution initially supports macOS/Linux.
- Added schema-1 project presets and a built-in `5w6w` planning policy, with
  strict TOML validation and checksummed local resources. Added `preset check`,
  `campaign plan` and read-only `campaign status`: campaigns freeze dataset
  selection, configuration, input hashes and portable model/preset snapshots,
  retain blocked datasets with diagnostics, report duplicate observations and
  detect changed or missing inputs after relocation. These commands plan input
  and model selection; scientific preflight and execution remain separate
  commands. Existing stage and checkpoint schemas are unchanged.
- Published the DOHU input dataset (configuration, authoritative MTZ, autoPROC/STARANISO
  CIF, and summary) for collaborator reruns. Added `examples/README.md` and an
  eight-file SHA-256 manifest covering DOHU and the existing dated QiC example.
  The uploaded short-name QiC inputs match that existing example byte for byte;
  no duplicate dataset or new run history is added. Input bytes are preserved.
- Added local dictionary construction for modified-residue mutations with one
  unambiguous NARestraints record and a supported DNA/RNA construction parent.
  Included `OHU.cif` supports the D:OHU example. Curated overrides take priority;
  missing dictionaries are not downloaded automatically.
- Added Phenix 2.2.x to the file-scoped Data Manager selector policy, alongside
  2.1.x; 1.20.x retains legacy explicit selectors. Other version families still
  stop before allocating refinement state.
- Added a two-stage PostMR phosphate audit. Verified internal OP3/O3P leaving
  atoms are removed before NARestraints and after ReadySet, authorized terminal
  oxygens remain (see explicit recipe/dataset policy above), raw products are preserved, and measured geometry/removals are
  recorded under `phosphate_cleanup`. Ambiguous or damaged phosphates stop
  preparation instead of receiving guessed coordinate repairs.
- AutoSol reports now record an optional checksummed density-map reference,
  availability status, and diagnostic independently of phase acceptance.

- Added a tracked POSIX `./nasolve` source-checkout launcher. It resolves
  symbolic links, uses the checkout's `.venv` by default, accepts an explicit
  `NASOLVE_PYTHON`, and scopes `src` path injection to the launched process.
- Added portable run, dataset, and repository artifact anchors. Newly created
  checkpoint registries no longer depend on a collaborator's checkout path.
- Added compatibility rebasing for legacy run reports containing absolute
  paths from another computer.
- Added `nasolve workspace use`, `status`, and `clear`. An active dataset/run is
  stored in the existing user configuration, allowing stage commands to omit
  repeated paths.
- Added collaboration and AI-review guidance in `AGENTS.md`,
  `docs/collaboration.md`, and the pull-request template.
- Added a frozen SHA-256 for authoritative reflections to new AutoMR reports.
- Added a checksummed, run-local `Model/seq_base.txt` snapshot for new standard
  frame runs when the selected catalogue supplies it, so conditional AutoSol
  can continue without the source catalogue.
- Added SHA-256 values for AutoSol heavy-atom and refinement-phase outputs.

### Fixed

- Checkpoint listing now validates and displays an uninitialized PostMR root
  without creating a registry or directories; existing registries are read
  without migration. Modifying checkpoint operations still initialize normally.
- Standalone AutoRefine prints a shell-quoted inspection command for its exact
  result checkpoint when the existing `show` resolver can resolve its model and
  maps, including unselected results and recorded map fallbacks. Resolution
  failures are reported explicitly, and failed results no
  longer receive a generic suggestion to select them.

- ReadySet exit status 0 with the explicit `No unknown residues` result is now accepted as a successful no-op when Phenix writes no `*.updated.pdb`. NASolve preserves the already-prepared model, still runs its phosphate/atom-count audits, and records the no-op output mode; missing ReadySet output without that explicit success condition still fails closed.

- Refine Doctor's interactive `i` answer now opens the recommended checkpoint
  in Coot and returns to the same selection prompt. Answers take effect after
  Enter; only `y` selects a checkpoint. Declining, blank input, EOF or Ctrl-C
  preserves the current selection, and an inspection-launch error allows a
  retry. Printed inspection/selection commands now quote paths with spaces.
- AutoSol now uses fully qualified PHIL parameter paths, including nested
  phasing, model-building and general settings. This fixes Phenix 2.2 rejecting
  the ambiguous `data=` argument before reading inputs. Standalone and campaign
  AutoSol retain the same MR model, anomalous labels, wavelength and no-building
  policy. The external-tool fixture now rejects abbreviated parameter names.
- Campaign worker shutdown now recognizes macOS process groups containing only
  exited, unreaped workers when group probes return `EPERM`. Cleanup still
  requires verified inactivity, reaps its owned child and checks the group
  afterward; live or unidentified groups retain the existing error and retry
  safeguards.
- Curated residues already present in a model retain their dictionaries and
  component-identity records even when no mutation is needed, including
  nonstandard modified-pair preparation.
- PostMR phosphate cleanup now tolerates reused atom serials when full atom
  identities are unambiguous. Atom-associated records are removed by identity
  within their model/TER segment; unrelated atoms sharing a serial survive.
  Affected CONECT records still stop preparation if their endpoint serials
  cannot be resolved uniquely. Phosphate geometry guards are unchanged.

- Corrected view/AutoSol test expectations for canonical paths on macOS, where
  `/var` aliases `/private/var`. The affected tests now use an explicit symlink
  fixture on every host, preserving relocation and checksum assertions while
  checking resolved output paths. Runtime path resolution is unchanged.

- `show` now follows the selected current checkpoint, including manual and
  PostMR nodes. Newer unselected attempts do not replace that view. Manual
  models inherit clearly labelled ancestor maps only with matching observations.
- PostMR and AutoSol views now load the prepared model and ligand dictionaries.
  PostMR prefers accepted AutoSol density when available; AutoSol displays its
  density-modified map and separate heavy-atom overlay. Console and launch
  records identify the model and map sources, with missing/corrupt declared
  sources reported explicitly.

- AutoRefine now omits unsupported Phenix Data Manager Miller-array blocks on
  Phenix 1.20.x while retaining the exact observation, Free-R, and optional HL
  phase file/label selectors on the command line. Phenix 2.1.x keeps the
  file-scoped Data Manager block that prevents ambiguous MTZ-array selection.
- Refine Doctor now passes the same version-derived reflection-selector policy
  to every bounded AutoRefine branch.
- Fresh-terminal source launches no longer depend on macOS Python loading the
  editable-install `.pth` file; this avoids `ModuleNotFoundError: nasolve` when
  that metadata carries the hidden file flag. The launcher also places the
  checkout package ahead of the current directory, so an unrelated
  `nasolve.py` cannot shadow it.
- Bound AutoRefine observation, Free-R, and experimental-phase Miller arrays
  explicitly to their source MTZ files through the Phenix Data Manager. This
  fixes Phenix 2.1 preflight failures when an MTZ contains multiple equally
  suitable observation arrays, such as `IMEAN,SIGIMEAN` and anomalous
  intensities.
- Committed PostMR runs can now resolve their prepared model, authoritative
  observations, restraints, Phaser outputs, and optional AutoSol inputs after
  the repository is cloned or moved.
- Legacy paths now prefer the selected checkout even when the collaborator's
  old checkout still exists, understand Windows/POSIX path forms, and reject
  traversal or ambiguous rebases.
- AutoRefine recovers approved AutoSol phases after relocation when the root
  checkpoint was created before AutoSol completed.
- AutoRefine now fails closed when an approved phase source is unavailable,
  instead of silently continuing without the requested phases.
- Malformed saved workspace, run-report, and checkpoint JSON now produces a
  guarded NASolve error instead of an uncaught type exception.
- Coot views now fail closed when a reported map, primary refinement MTZ, or
  required ligand dictionary is missing or fails checksum validation.
- AutoSol and Coot now verify available frozen model/data checksums before
  consuming relocated artifacts.
- Legacy AutoSol runs recover the exact declared frame's `seq_base.txt` after
  relocation even when the historical pair-model filename no longer exists;
  the selected checkout wins while an older checkout is still present.
- Coot scratch creation rejects unsafe checkpoint path components and existing
  symlinks that could redirect logs, launch records, or backups outside a run.

### Changed

- Every AutoRefine checkpoint and round report now records the discovered
  Phenix version and the selected `legacy-explicit` or
  `data-manager-file-scoped` reflection-selector mode. Refine Doctor reports
  record the same provenance. Unknown, malformed, and unvalidated Phenix
  version families stop before allocating a refinement or Doctor directory.
- AutoRefine reuses the checksum produced by mandatory phase validation when
  writing the child checkpoint, avoiding a second full read of the AutoSol MTZ
  in each round. A stable file-identity check before and after Phenix rejects
  phase replacement, modification, or deletion during the external process.
- AutoRefine checkpoint serialization hashes each unique output file once, so
  legacy Phenix layouts that use one MTZ for maps and refinement reflections
  no longer read that potentially large file twice.
- AutoRefine now writes file-scoped Miller-array selections into
  `autorefine.params` while retaining the existing `xray_data` refinement
  parameters for compatibility.
- Separate AutoSol phase files receive their own file-scoped
  Hendrickson-Lattman label selection.
- `automr` uses the active workspace dataset when its dataset argument is
  omitted. `postmr`, `autosol`, `autorefine`, `refine-doctor`,
  `checkpoints list`, `checkpoints add`, `checkpoints use CHECKPOINT`, and
  `show` can use the active run.
- Checkpoint schema 2 uses `anchor` plus `relative_path` (or an explicit
  `absolute_path` fallback), SHA-256, and size. Refinement review outputs such
  as maps, reflection/mmCIF files, metrics, and logs use the same references.
- A newly allocated AutoMR run becomes active automatically. Explicit one-off
  paths passed to later stages do not silently change the active workspace.
- Machine-local Coot GUI scratch and PostMR Coot backup directories are ignored;
  portable stage artifacts remain visible for deliberate commits.

### Tests

- Doctor triage regression coverage includes the QiC inversion ranking case,
  ordinary IMEAN fallbacks, fixed and f''-only scattering, alternate selections,
  missing prerequisites, unavailable audits, failed external runs, budgets,
  immutable sibling lineage and inspection without selection. The complete
  local suite passes 393 tests and 129 subtests, with two existing platform
  skips. New Phenix recipes still require live validation on the user runtime.
- The assembled fixes pass 225 tests. The full local run with
  `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src` also passed 51 subtests; the user
  confirmed 225 passing tests on the Apple Silicon checkout.
- On 2026-09-09, DOHU `AutoMR/run_004` completed AutoMR, PostMR, and a full
  five-cycle AutoRefine execution with Phenix 2.2.1-6174 and Coot 1.3.3:
  TFZ 14.30, final Rwork/Rfree 0.147/0.209, and current checkpoint `refine-001`.
  The user confirmed successful Coot inspection. The
  [validation record](docs/validation-dohu.md) preserves the raw clashscore,
  lattice context, execution details, and limits of this check.
- Added coverage for current/manual checkpoint views, AutoSol map provenance,
  relocated and symlinked paths, and phosphate cleanup with duplicate serials
  and atom-associated records.
- Added Phenix 1.20.1-4487 compatibility regressions that emulate rejection of
  the Data Manager block and verify exact anomalous/mean observations, Free-R
  flags, separate HL phases, automatic target selection, Doctor propagation,
  and report/checkpoint provenance without production use of `--unused_ok`.
- Added fail-before-mutation coverage for unknown or unvalidated Phenix
  versions and retained the Phenix 2.1 file-to-array binding regressions.
- Revalidated the post-patch `data-manager-file-scoped` path from a relocated
  temporary Q:iC run with Phenix 2.1-6048. Dry-run preflight and a one-cycle
  refinement both exited successfully; the intentionally short refinement was
  retained as `AUTOREFINE_REVIEW` at Rwork/Rfree 0.177/0.160.
- Added source-launcher coverage for fresh terminals, symlink resolution,
  argument and working-directory preservation, `PYTHONPATH` scoping, explicit
  and default interpreters, current-directory import shadowing, and
  missing, directory-valued, and non-executable interpreter diagnostics.
- Added hash-call regression coverage for inherited/recovered AutoSol phases
  and shared AutoRefine MTZ outputs, plus phase drift during Phenix, immutable
  failed-attempt retention, and successful next-round recovery.
- Added relocation coverage for reports created beneath a collaborator's
  absolute directory.
- Added coverage for anchored run/dataset references and path-escape rejection.
- Added repository-anchor and exact frame-sequence relocation coverage,
  including old and current checkouts existing at the same time.
- Added same-machine-copy, cross-platform legacy-path, schema migration, and
  checksummed-reflection coverage.
- Added legacy AutoSol coverage for reports that predate frozen input hashes.
- Added workspace configuration and optional-target parser coverage.
- Added regression coverage for mean-intensity refinement without AutoSol.
- Expanded coverage for anomalous observations, Free-R flags, and a separate
  experimental-phase MTZ.
- Validated the AutoRefine selector patch with Phenix 2.1-6048 in both a
  `--dry_run` preflight and a complete five-cycle E:G refinement that reached
  `AUTOREFINE_READY`.
- Re-ran all four example workflows end to end with Phenix 2.1-6048 and Coot
  1.3.1. Fresh five-cycle refinements reached `AUTOREFINE_READY` for E:G
  (TFZ 11.80, Rwork/Rfree 0.127/0.168), the F:D synthetic integration fixture
  (TFZ 11.80, 0.131/0.171), Q:E (TFZ 10.80, 0.161/0.208), and Q:iC (TFZ
  10.80, 0.146/0.148). Q:iC retained the expected non-blocking AutoSol warning
  and refined its anomalous iodine without experimental phases.

### Migration

- New phosphate checks apply to newly prepared PostMR outputs; existing models
  and immutable runs are not rewritten. Use a new run to reproduce preparation
  with these checks. The new AutoSol map reference is optional for older reports.
- No configuration or artifact migration is required for the source launcher;
  the generated package entry point remains supported when its editable
  installation metadata loads normally.
- Existing absolute-path reports remain readable; they do not need to be
  rewritten before use.
- Existing checkpoint schema-1 registries are read in place and migrate to
  schema 2 on their next write. Older NASolve versions reject schema 2 instead
  of interpreting an anchor-relative path against the current directory.
- Portably resolvable legacy references are normalized during migration;
  unavailable or external-only legacy evidence is retained explicitly as
  non-operational provenance rather than remaining a raw schema-1 path in a
  schema-2 registry.
- Legacy AutoMR reports without `reflections_sha256` remain usable and are
  marked `legacy-unverified` when a root checkpoint is first initialized,
  rather than silently claiming observation checksum verification.
- Existing failed refinement rounds remain immutable. Running
  `nasolve autorefine RUN` again creates the next numbered checkpoint.
- Workspace state is machine-local and is not intended for Git.

### Known compatibility notes

- Phenix 1.20.x uses explicit file/label command selectors without the newer
  Data Manager block; Phenix 2.1.x and 2.2.x use both. User validation also
  completed Q:iC refinement and bounded Refine Doctor branches on 1.20.1-4487
  with `legacy-explicit` recorded in the refinement report.
- File-scoped Data Manager parameters remain validated with Phenix 2.1-6048
  in preflight and complete refinement runs.
- The new live Phenix 2.2.1-6174 check covers mean-intensity refinement with
  no experimental phases and anomalous refinement off. It does not constitute
  a new live AutoSol or anomalous-refinement validation for that version.
- Large MTZ/map histories may eventually warrant Git LFS; current fixtures are
  intentionally kept in ordinary Git until a repository-wide LFS policy is
  adopted.

## [0.1.0 development baseline] - 2026-08-28

### Added

- External Phenix and Coot discovery and saved runtime configuration.
- Guarded AutoMR with frozen inputs and TFZ classification.
- PostMR model preparation, curated modified nucleotides, NARestraints, and
  hydrogen-free ReadySet.
- Conditional AutoSol for supported nucleotide-bound I, Br, and Se atoms.
- Sequence mutation, modified-pair restraints, and model mirroring.
- Checkpointed AutoRefine and bounded Refine Doctor triage.
