# Changelog

All notable user-visible, compatibility, schema, and reproducibility changes
to NASolve are recorded here.

The format follows Keep a Changelog. Until the first tagged release, historical
entries are reconstructed from repository history.

## [Unreleased]

### Added

- `nasolve autorefine --geometry [METAL=]PRESET` (user's specification,
  2026-11-24): a preset coordination geometry for a `--coordination` site.
  For now `sqp`, square planar: four donors, given in order around the
  square; one angle restraint donor-metal-donor (metal at the vertex) for
  each pair given next to each other, the last with the first, 90 deg, sigma
  6 (`nasolve.metal_coordination.GEOMETRY_PRESETS`). Pairs across the square
  are not restrained. A site with another number of donors is refused; donors
  given next to each other but more than 135 deg apart in the starting model
  are reported (probably not in order). A bare preset applies to the only
  `--coordination`; with several, `METAL=PRESET` names the site. Without
  `--geometry` nothing changes (linear angle for two donors, none for more).
  Site records gain `preset` and `angles` (all angle restraints; `angle`
  still holds the linear one). phenix.refine simulated.
- `nasolve autorefine RUN --coordination METAL=DONOR,DONOR` (user's
  specification, 2026-11-23): restrains a heavy-atom site to its binding
  partners with phenix.refine `geometry_restraints.edits`
  (`nasolve.metal_coordination`). One bond per donor with the lab's lengths
  for now (Ag-O 2.5, Ag-N 2.3, Hg-O 2.1, Hg-N 2.1 A), sigma 0.2 A; with two
  donors (linear) also the angle donor-donor-metal (vertex at the first
  donor) 0 degrees, sigma 3, a non-physical angle that keeps the metal on the
  donor-donor line; no angle for three or more donors. Written to
  `round_NNN/metal_coordination.eff`, passed to phenix.refine, recorded in the
  round report (`inputs.metal_coordination`: sites, current distances and
  angle, warnings) and the checkpoint (`metal_coordination`), and kept by
  later rounds; restraining a metal the parent checkpoint already restrains
  is refused before a round exists. Atoms are checked in the round's starting
  model (first model; alternate conformations and quoted atom names refused
  for now). The file format was parsed against Phenix 2.2.1's own edits
  definition; cctbx handles the 0-degree angle (its gradient guard for
  collinear atoms). Phenix's `phenix.metal_coordination` writes the same
  format but its lengths (Harding 2006) do not cover Ag or Hg. phenix.refine
  simulated.
- `nasolve autofdp --nudge Q` (user's request, 2026-11-19): the metal site
  starts the f'' rounds at its occupancy in the model (AutoRefine's) plus Q
  (two decimals; negative lowers it; must stay within `--occupancy-floor` and
  1, else refused before anything is written). One site for now
  (`autofdp.nudge_targets`, to become a selection later): the model's only
  metal site, or the one `--site` names. `start_model.pdb` gets the occupancy
  (and B) change; `report.json` records `settings.nudge` and each site's
  `occupancy_start` and `nudge`; the plan and `autofdp.log` show
  "occupancy 0.90 -> 0.95 (--nudge 0.05)". First written as minus Q; the
  user corrected it to plus Q the same day. phenix.refine simulated.
- `nasolve postmr RUN --merged-model FILE`: start PostMR from a chosen model
  with merged metal sites instead of the newest merge (user's request,
  2026-11-19). FILE is a name in the run's `Phaser/` folder (e.g.
  `mr_solution-merged-2.pdb`) or a path. A recorded, unchanged merge keeps its
  record. Any other PDB file is used as it is: a model saved by hand in Coot,
  or a recorded name whose file changed since (U-Y_Ag run_004 has both: the
  first merge renamed `mr_solution-mergedadad.pdb`, and a two-Ag model saved
  as `mr_solution-merged.pdb`, which the default refuses as changed). Such a
  file is recorded with `recorded_merge: false`, its checksum and
  `comparison_with_phaser_model`, residue by residue: residues added (the
  merged metals), residue-name changes (marked planned when they are the run's
  planned mutations, which the metal search builds into its partial model),
  and atoms added, removed or moved (> 0.01 A) in the other residues. A warning
  is printed and logged for it, and another when residues were removed or
  changed outside the plan, or atoms removed or moved. On the real U-Y_Ag
  run_004 `mr_solution-merged_2Ag.pdb`: Ag B:101 and B:102 added, A:12 DG->DU
  (planned), nothing else. A missing file or one with no atoms is refused before anything is
  written (the recorded merges are listed). Excludes `--original-model`. The
  run's own selection is not changed; `metal_merge.selected_by` records how
  the model was chosen, and the command prints the other recorded merges.
- `nasolve autofdp ROUND`: AutoFDP from an AutoRefine round folder
  (`RUN/AutoRefine/round_NNN`; user's request, 2026-11-19). Its
  `refined_NNN.pdb` is the starting model, the one X-ray data file its
  `refined_NNN.eff` names (`data_manager.miller_array` amplitude or intensity
  arrays; Free-R flags and phases are not data) is the MTZ, and that .eff is
  the final round's (as `--final-eff`; giving both is refused). The round's
  `report.json` must name the same model; its checkpoint and status are
  printed and recorded in the run report as `inputs.autorefine_round`
  (anchored references with checksums of the model and .eff, and the .eff's
  data labels). When the f'' rounds' anomalous array differs from the .eff's
  data labels, a warning is printed, logged and recorded (`warnings`). The
  final-round edits were checked on the real Y-T nasolve run_002
  `round_002/refined_001.eff`: its model, Ag group and `group_anomalous` are
  removed and the round's `autorefine.params` is left out, the restraint
  PHIL files are kept. The two-file form now refuses extra arguments
  (a third path was ignored). phenix.refine, ReadySet and Coot simulated.
- `nasolve postmr RUN --force-overwrite`: rerun PostMR on an AutoMR run that
  already has one. Since 2026-11-23 (user's request) it keeps the previous
  result instead of deleting it (it deleted from 2026-11-19): after every
  check has passed, `PostMR/`, `AutoSol/`, `AutoRefine/` (rounds and
  checkpoints), `RefineDoctor/` and `CootGUI/{postmr,autosol,autorefine}` are
  moved unchanged to `PostMRHistory/postmr-NNN/` (`nasolve.postmr_history`),
  their run-report records saved in its `result.json`, the report reset to
  the Phaser result and PostMR run again. Each PostMR result is numbered
  (`postmr_result.id`, `postmr-001` onwards). `nasolve postmr RUN --history`
  lists the results (time, starting model, AutoRefine rounds and current
  checkpoint with R values, stage and status); `--switch postmr-NNN` keeps
  the current result and moves the chosen one back with its records
  (`postmr_switches`), so every path it recorded is valid again. Phaser
  outputs, metal searches and merges belong to the run. The Phaser/TFZ gate
  is applied again (so a refused rerun moves nothing).
  Schema: the run report gains `phaser_result` (the Phaser status and message
  PostMR replaces with POSTMR_READY; for older reports re-derived from
  `execution.phaser`) and `postmr_overwrites` (each rerun: the number the
  previous result was kept as, the folders moved, the previous stage, status
  and PostMR time). The AutoMR status logic is now
  `phaser.classify_phaser_result` / `classify_phaser_record` (same rules and
  messages). A rerun that fails leaves the run at its Phaser result. ReadySet
  simulated in tests.

- PostMR runs NARestraints on runs without a frame (nonstandard mode): its
  guesser finds the model's base pairs and PostMR restrains all of them, plus
  stacking, with NARestraints' defaults (canonical pairs; metal-mediated pairs
  detected and not restrained, nor their bases' stacking). Previously such a
  run got no pair or stacking restraints (`narestraints.mode = "none"`), e.g.
  Y-T nasolve AutoMR run_001. New files `PostMR/Restraints/guessed_pairs.txt`
  and `narestraints_guessed_pairs.phil` (passed to AutoRefine); report
  `narestraints.mode = "guessed-pairs"` with candidates, warnings and
  `metal_mediated_pairs`. `nasolve postmr --no-narestraints` keeps the old
  behaviour; it is refused with a frame or `--modified-pairs-only`. The
  terminal and `postmr.log` gain a one-line NARestraints summary for every
  mode. On the real `MR_frames/5W6W/5W6W_noPO4.pdb` (NARestraints d1bb52d):
  19 pairs, 38 stacking pairs. ReadySet and phenix.refine simulated.

- `nasolve autorefine --restraints FILE` (repeatable): extra restraints for a
  round, PHIL (`.eff`, `.phil`, `.params`) or CIF (`.cif`); other file types
  are refused. Each file is copied into `round_NNN/extra_restraints/`, passed
  to phenix.refine as that copy and recorded in the checkpoint's restraints
  (run-anchored, SHA-256), so later rounds from that checkpoint keep it; a
  file whose contents are already in the parent checkpoint's restraints is
  refused before a round exists. The round report's `inputs.user_restraints`
  lists each original path, checksum and copy; the terminal shows them.
  Motivation: Y-T nasolve AutoMR run_001 (nonstandard, no frame) had no
  secondary-structure restraints, and a `5W6W_secondary_structure.eff` copied
  into `PostMR/Restraints` afterwards was not used. phenix.refine simulated.

### Changed

- The anomalous map opens in Coot as an ordinary map, named "Anomalous
  (ANOM/PANOM)", not as a difference map (user's request, 2026-11-23): gold
  (`set_map_colour` 1.0, 0.84, 0.0; `coot_view.ANOMALOUS_MAP_COLOUR`),
  positive contour at 3 sigma. An anomalous map already open (from
  the user's Coot settings) is switched the same way. All NASolve Coot views.
- AutoFDP accepts a refined f'' within +/-10% of the theoretical f'' of the
  site's element by default (user's request, 2026-11-23), instead of a fixed
  +/-1 e- (`--tolerance`, removed). `--acceptance_absolute E` (also spelled
  `--acceptance-absolute`) gives a fixed window of E electrons for every
  site. Each element now has its own window: Ag (2.03 e- at 1.0074 A)
  +/-0.20 e-, Hg (10.1 e-) +/-1.01 e-. With +/-1 e-, an Ag site (f'' about
  2 e-) was accepted at any occupancy from about 0.5 up, so runs that started
  at 1.00 ended at 1.00 (7CUA: the same data were accepted at 0.75 and at
  1.00); Hg runs are judged about as before. The plan, the command's tables
  ("within +/-0.20", per site in e-), `autofdp.log` ("Acceptance: refined
  f'' within +/-10% of theory (AG +/-0.20 e-)") and the run reasons ("within
  +/-10% (0.20 e-) of theory") show the window. **Schema:** AutoFDP
  `report.json` schema 3: `settings.acceptance` (`rule` `relative` with
  `fraction` 0.1, or `absolute` with `absolute_electrons`) replaces
  `settings.tolerance`, and each round site's `acceptance` (e-) and
  `within_acceptance` replace `within_tolerance`; earlier runs are not
  rewritten. **Compatibility:** scripts passing `--tolerance E` must use
  `--acceptance_absolute E`; `fdp_round` and `execute_autofdp` take
  `acceptance_absolute` instead of `tolerance`. phenix.refine simulated.
- PostMR keeps the bases of metal-mediated pairs stacked on their chain
  neighbours (user's request, 2026-11-19): NARestraints still writes no
  base-pair restraints for such a pair, but each base keeps its stacking
  restraints (`keep_metal_pair_stacking=True` passed in frame-template and
  guessed-pairs mode when NARestraints has the option; recorded as
  `narestraints.metal_pair_stacking = "kept"`; the PostMR summary says
  "(stacking kept)"). Before, NARestraints left those bases unstacked and
  stacked nothing across them: in Y-T nasolve run_002 (Ag between A:12 DYI
  N2T and B:4 DT O4) A:11-A:12, A:12-A:13, B:3-B:4 and B:4-B:5 had no
  stacking; with this change the real PostMR model gets all four (A:12's as
  parallelity edits). NARestraints fork `feature/anomalous-scatterers`
  b6993e9 makes the same the default there. Runs already prepared keep their
  restraints; rerun PostMR (`--force-overwrite`) to get it.
- DYI dictionary: the central C4-N4 bond is free to rotate (user,
  2026-11-23): each ring has its own plane again (`plan-1`, dC's base plane
  with N4; `plan-2`, the triazole with C4), and TZF's ring-to-ring torsion
  (C5-C4-N4-C5T, 180 deg, period 2, e.s.d. 5) is left out, so nothing
  restrains rotation about C4-N4. The curated check enforces it
  (`CuratedLigand.free_rotation`: no torsion about the bond, no plane with
  atoms on both sides). The ideal coordinates still have the triazole
  coplanar (one conformation; moved by at most 0.001 A). Reverses the merged
  plane below. Runs prepared earlier keep the dictionary PostMR recorded;
  rerun PostMR (`--force-overwrite`) to use this one, and do not pass
  `Y-T nasolve/DYI_A12_coplanar.eff` any more.
- DYI dictionary (superseded on 2026-11-23, see above): one base plane over
  both rings (user's request, 2026-11-19). dC's base plane and the triazole plane are merged into `plan-1`
  (17 atoms: both rings, C1', O2, the ring hydrogens; e.s.d. 0.020 A), so
  refinement keeps the triazole in the plane of the pyrimidinone. With a plane
  per ring, sharing only C4 and N4, the triazole turned 38 deg out of the base
  plane in Y-T nasolve AutoMR run_002 refine-001 (the ring-to-ring torsion,
  e.s.d. 5 deg, was the largest dihedral residual). Built by
  `scripts/build_dyi_dictionary.py` from the same Monomer Library files; ideal
  coordinates move by at most 0.001 A. The curated check now requires one
  plane holding both rings. Runs prepared before this keep the dictionary
  PostMR recorded; rerun PostMR (`--force-overwrite`) or add a planarity edit
  with `autorefine --restraints` to use it. Not yet run in phenix.refine.
- Metal searches (`nasolve autosol RUN --metal ELEMENT`) are written inside the
  AutoMR run they come from, `AutoMR/run_NNN/AutoSolMetal/run_NNN` (numbered
  per AutoMR run), instead of `DATASET/AutoSolMetal/run_NNN` (user's request,
  2026-11-18). Searches already in `DATASET/AutoSolMetal` are not moved and
  still open, from the folder or with `--stage metal`. Schema: new index
  entries in the AutoMR report (`metal_autosol.searches[].directory`,
  `.report`) and the merge record's `metal_search` are run-anchored
  (`AutoSolMetal/run_NNN`); the search's own report keeps dataset-anchored
  references to the AutoMR run and inputs. `run_context.dataset_directory`
  recognises `AutoMR/run_NNN/AutoSolMetal/run_NNN`. AutoSol simulated in tests.

### Fixed

- The anomalous difference map opens in Coot again (reported 2026-11-23:
  after AutoRefine rounds only the 2mFo-DFc and mFo-DFc maps opened). Two
  causes in the startup script: Coot (bundled with Phenix 2.2.1) requires a
  Python bool for `valid_labels`' use-weights argument and raised "argument 5
  of type 'bool'" on the 0 it was given, so the script stopped at once
  (`NASOLVE_ANOMALOUS_MAP ... "status": "error"` in every U-Y_Ag run_004
  `CootGUI/autorefine/*/coot_gui.log`); and it looked only for `ANOM`/`PHANOM`
  (phenix.maps' labels), while phenix.refine writes `ANOM`/`PANOM`. It now
  passes `False` and tries `ANOM`/`PANOM`, then `ANOM`/`PHANOM`; the map is
  named after the pair found and the log line records it. `launch.json`
  `anomalous_map.labels` is now the list of pairs tried. Applies to every
  NASolve Coot view (show, AutoRefine, Refine Doctor's `i`, AutoFDP final
  round). Checked against Phenix 2.2.1's `coot.py` signatures and the real
  MTZ columns; not run in a live Coot here.
- PostMR and AutoFDP give phenix.ready_set one dictionary file per curated
  residue (`ready_set.input.cif_file_name` repeated;
  `ligand_profiles.readyset_dictionary_files` splits a file defining several).
  ReadySet (Phenix 2.2.1, `ready_set_utils.get_user_defined_restraints`) reads
  only the first component of each file, so with NASolve's one combined file
  every other curated residue counted as unknown. In U-Y_Ag nasolve AutoMR
  run_004 (DU and DYI, 2026-11-19) ReadySet used the file for DU, fetched the
  wwPDB CCD entry DYI, which since 2026 is an unrelated molybdopterin
  ("Molybdopterin hydrox dioxo molybdenum", C10 H13 Mo N5 O9 P S2; not in
  Phenix 2.2.1's bundled CCD), ran eLBOW on it with the model's atoms and
  stopped ("Bond length too long : 10.47"). The same cause gave DE its
  from-scratch restraints beside 1AP. The combined `curated_ligands.cif` /
  `input_dictionaries.cif` is still written and checked. ReadySet simulated
  (a stand-in that reads dictionaries as Phenix 2.2.1 does).

- PostMR stopped with "ReadySet failed with status 1" (eLBOW:
  `AttributeError: 'AtomClass' object has no attribute 'atom_name_in_input'`,
  Phenix 2.2.1) when a curated residue already in the starting model used its
  own atom names: Y-T_Ag AutoMR run_005, whose search model `DYI_G.pdb` has a
  DYI with C1, C3, N2, N5, O1, ... instead of DYI.cif's C5, C2, N3, N4T, O2,
  .... PostMR now checks every curated residue it keeps from the starting
  model before `PostMR/` exists: each heavy-atom name must be in the
  dictionary, and each dictionary bond between two of its atoms must be a bond
  in the model (so names on the wrong atoms are caught too). The error names
  the atoms and the dictionary's names. A dictionary without an atom list gets
  only the bond check. Tested on the real `DYI_G.pdb` (refused) and on a
  renamed copy (accepted); ReadySet itself not rerun.

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

- Curated ligand `DYI`: O4-triazolyl-dU (the suppliers' name for
  4-(1,2,4-triazol-1-yl)pyrimidin-2-one 2'-deoxyriboside, a convertible
  nucleoside), under a laboratory label; it is not a CCD component (checked
  against the CCD in Phenix 2.2.1-6174, entries to 2026-07-10). Use the code
  itself (`--pair DYI:G`, `A:12 = DYI`); there is no alias. Parent `DC`;
  the triazole's N1 is `N4`, the other ring atoms `N2T`, `C3T`, `N4T`, `C5T`.
  NARestraints sees it as `C34` (N4-methyl-dC). `DYI.cif` is a reviewed
  restraint dictionary built by the new `scripts/build_dyi_dictionary.py` from
  the CCP4 Monomer Library AceDRG entries `DC` (shared atoms, unchanged) and
  `TZF` (triazole ring), commit 0af022b, with C4-N4 1.411 +/- 0.015 A; its
  coordinates satisfy its restraints (bonds within 0.011 A, angles 1.5 deg).
  PostMR places `N4` on the parent's C4-N4 line at that length. Like `1AP`,
  `DYI` is used for refinement in place of ReadySet's dictionary (ReadySet
  reads one component from NASolve's combined dictionary file and builds any
  other non-CCD residue from the model, without chirality restraints, as it
  did for DE beside 1AP). Not yet run in Phenix or Coot (simulated in tests).
- Reviewed restraint dictionaries are now declared per curated entry
  (`reviewed_restraints`, `required_planes`); `AUTHORITATIVE_CODES` is derived
  from them, and the 1AP-only checks apply to every such entry, now also
  requiring signed sugar chiral centres. Schema: the PostMR precedence policy
  is `reviewed-dictionaries-over-readyset-v1` (was
  `reviewed-1ap-over-readyset-v1`; nothing reads it, and old runs keep theirs).
- `--tfz-override` for `nasolve postmr` and `nasolve autosol RUN --metal
  ELEMENT`: continue from a Phaser result whose only problem is a TFZ below
  8.0 (`MR_REVIEW`, or `MR_FAILED` below the 7.0 floor). Phaser must have
  exited with status 0, reported a TFZ, written both the solution model and
  map, and used the requested reflection data; any other failure is still
  refused, with its reason. The 7.0/8.0 cut-offs and the run's MR status do
  not change, and campaigns do not use the option. Schema: the PostMR record
  (`postmr.tfz_override`, `PostMR/report.json`) and the metal-search report
  (and its entry under `metal_autosol.searches`) gain `tfz_override`: null,
  or the accepted TFZ, original status and message, the cut-offs and the time.
  The command prints a warning, and `postmr.log` records it. Simulated in
  tests; not yet run on real data.
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
    not move. `--theory ELEMENT=FDP` replaces an element's theoretical f''.
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
    still on the wrong side. The last round's model is `final_model.pdb`;
    the command prints every adjustment and a final table. Status `AUTOFDP_CONVERGED` (exit 0) when every calibrated
    site is within tolerance, `AUTOFDP_REVIEW` (exit 2) with per-site reasons
    otherwise, `AUTOFDP_FAILED` (exit 1) when a round fails. Schema (additive,
    still 2): each `rounds` entry gains `input_model`, `adjusted_site`,
    `occupancy_before`/`after`, `method` and `note`, and each site its
    `alternate`; `report.json` gains `occupancy` (method, limits, every
    adjustment, per-site outcome and history), `final_model`, and
    `settings.max_iterations`/`occupancy_floor`; `steps.occupancy_adjustment`
    is now `done`, `review`, `failed` or `not run`. The `decide_site_update`
    stub is replaced by `propose_occupancy`.
  - Step 3, f'' moving with occupancy: when an adjustment lowers a site's
    occupancy and its f'' falls (or raises it and f'' rises), the site goes
    back to its earlier occupancy with the highest (lowest) f'' in a restore
    round and is held there while the other sites are adjusted; it is taken
    up again when no other site can be adjusted (same per-site limit), and a
    second such response stops it there for review. Before, the proportional
    step kept lowering it (PdUHgPdU Hg h:1, live on Phenix 2.2.1: f'' peaked
    near occupancy 0.3 and the site was pushed to the 0.05 floor while its
    neighbours' f'' rose to 12-15 e-). New `wrong_response`; rounds may have
    `method: "restore"` (its note says why); outcomes gain `held`; restores
    are not counted as adjustments.
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
    `nasolve.autofdp`); `--element` is refused. phenix.refine, phenix.mtz.dump
    and phenix.python are simulated in tests; first run on real data
    2026-11-13 (PdUHgPdU, Phenix 2.2.1).
  - AutoFDP no longer calls NARestraints: it finds metal sites and their
    neighbouring nucleotides with the new `nasolve.metal_sites`, which holds
    NARestraints' heavy-atom search (`find_heavy_atoms`) and metal-pair search
    (`read_atoms`, `find_metal_mediated_pairs`) copied unchanged from
    AcetylCholine233/NARestraints `main` d1bb52d, plus `scan_metal_sites`
    (PostMR's selection: NASolve's metal list without Li/Be/Na/Mg/Al/K/Ca,
    outside nucleotides). Results are identical; tests compare them with
    NARestraints and PostMR whenever NARestraints is installed, and run
    AutoFDP with NARestraints unimportable. The neighbour method label is
    now "bridged base pair (metal-pair search)". PostMR and AutoRefine still
    use NARestraints as before.
  - No Coot for the f'' rounds: AutoFDP no longer opens an f'' round's model
    and maps in Coot (it did up to a70def8 on `ginkgo`), because those maps
    are not diagnostic. Only the final round (below) opens Coot.
  - Final round (`--final-eff MANUAL.eff`; partially implemented: it also
    serves as the B confirmation, and how its starting parameters are
    given may change; first live run 2026-11-18, PdUHgPdU AutoFDP run_013,
    Phenix 2.2.1, R-work 0.171, R-free 0.260): when every calibrated site
    converges, a phenix.refine round in `run_NNN/final_round/` with every
    setting of a manual refinement's .eff (for example the GUI's
    `.phenix/project_data/refine_NN.eff`), with only these changes: AutoFDP's
    `final_model.pdb` instead of the .eff's model (a reference model stays);
    the .eff's anomalous groups replaced by one group per metal site with f'
    at theory and f'' at its last refined value, `refine = None`, and
    `group_anomalous` unselected (Phenix 2.2.1 applies group f'/f'' without
    that strategy); every metal site added to
    `occupancies.remove_selection` (occupancy held; the .eff's own
    exclusions kept); outputs named `final_001`. The edited copy
    (`final_input.eff`, with a header) and the metal settings
    (`autofdp_final.params`) are written beside the outputs; the .eff
    itself is never changed. Extra PHIL files the .eff loads are checked:
    a refinement job's full parameter file is left out (phenix.refine's
    output .eff names the GUI's `refine_NN.eff`, which brought the old
    model back: "Wrong number of models of each type supplied", PdUHgPdU
    AutoFDP run_011), with its `default_phil`; another file that sets
    anomalous groups, a model, `group_anomalous` or the output name is
    replaced by an edited copy in `final_round/phil_files/`
    (`edits.removed_phil_files`, `removed_default_phil`,
    `copied_phil_files`; `FinalEff.job_phil_files`). Any other parameter
    naming the replaced model file is removed as well: phenix.refine loads
    every existing file a path parameter names, and the GUI's
    `refinement.gui.migration.refinement.input.pdb.file_name` brought the old
    model back as a second model (PdUHgPdU run_012;
    `edits.removed_model_references`). New line-preserving PHIL scanner
    `nasolve.phil_text` makes these edits without touching any other line.
    Preflight: the .eff must be a phenix.refine parameter file naming a
    reflection file, and every file it names must exist; another
    reflection file than AutoFDP's MTZ is a warning. Afterwards the held
    occupancies and the f'' shown in the log are checked. The result opens
    in Coot (`--no-coot`, `--coot` apply to it; `coot_view.launch_coot_files`
    is back). A failed final round makes the run `AUTOFDP_FAILED` (exit 1);
    a run that did not converge skips it. The `--final` stub (all
    parameters refined together) is replaced; `AutoFDPSettings.final` is now
    `final_eff`. Schema (additive, still 2): `final_round`,
    `settings.final_eff`, `steps.final_round`; new round status
    `AUTOFDP_FINAL_DONE`. Tested with simulated Phenix programs and one live run.
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
