# Data dictionary

- `05_data/internal16/`: topology/construct registry; experimental endpoint summaries; phenotype rank freeze; final PSCI values and exact dictionary; MD QC; 24 unordered pairs; six-dimensional system context.
- `05_data/literature_prior/construct_architecture_91.csv`: architecture covariates for 91 literature constructs.
- `phenotype_records_102.csv`: the 102 frozen literature phenotype records, including explicit missingness/evidence flags; no values were added.
- `validated_pairs_381.csv`: all validated literature pair labels before prior eligibility filtering.
- `eligible_non_tie_pairs_377.csv`: exact non-tie eligible prior set, with normalized direction and architecture summaries.
- `overlap_exclusion_audit.csv`: fold-specific audit governing direct held-out-system exclusion.
- `05_data/model/frozen_all_model_predictions.csv`: complete frozen V2 results; six model/sensitivity variants. Core reproduction uses M0–M3 only.
- `M2_FOLD_COEFFICIENTS.csv`, `M2_FOLD_SCALERS.csv`: exact fold artifacts reconstructed and verified in the provenance patch.
- `reported_model_comparison.csv`, `reported_system_metrics.csv`, `reported_construct_rankings.csv`: frozen reference summaries.
- `05_data/figures/`: exact source tables used/planned for computational manuscript panels. These retain negative and limiting evidence.
- Missing replicate-level experimental data remain missing. RHLARHLB_T04 has no numeric phenotype and is represented by its frozen rank/below-reference status only.

`D` follows the stored executable lineage: mass-weighted enzyme all-atom COM distance normalized by reference-coordinate enzyme CA radii of gyration. `L` follows the stored executable lineage: linker CA end-to-end distance divided by the same-frame sum of 16 consecutive linker CA distances. Older prose variants are retained as qualifications in the exact feature dictionary and were not used to alter values.

For manuscript-facing and executable definitions, `05_data/internal16/PSCI_EXACT_FEATURE_DICTIONARY.csv` and `03_code/psci/descriptors.py` are authoritative. `02_config/PSCI_PHENOTYPE_BLIND_FEATURE_FREEZE.json` is retained as historical phenotype-blind provenance and must not override the final executable definitions. Accordingly, final `D` uses mass-weighted enzyme all-atom COM distance normalized by reference CA Rg, and final `L` uses linker end-to-end distance normalized by instantaneous summed CA-CA contour length. The historical CA-only COM wording and fixed `(n_residues-1)*0.38 nm` contour approximation are not final executable definitions.

Missingness is preserved rather than filled: `F3` is not uniquely recoverable, `Q_topology` is not numerically instantiated, and `A`/`O` retain evidence-qualified missingness where applicable. No missing value is fabricated. The final predictive PSCI feature set is `[D, C, F1, F4, F5, L]`.


## ACS manuscript-facing additions

- `05_data/internal16/author_verified_reporting_annotations.csv`: author-verified `(G4S)3` experimental linker and `n = 3` reporting annotation; no individual replicates or SDs are invented.
- `05_data/figures/`: rebuilt to match the current ACS manuscript Figure 1–6 numbering and legends.
- `DATA_LICENSE.md`: software/data/literature-metadata licensing boundaries.
- `08_docs/PORTABILITY_AUDIT.md`: exact pinned-environment and compatible-environment numerical reproduction policy.
