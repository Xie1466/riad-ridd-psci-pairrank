# Methods-to-code map

| Method | Public code | Frozen data/config |
|---|---|---|
| PSCI extraction lineage and formula assembly | `03_code/psci/frozen_source/`; `03_code/psci/descriptors.py` | `PSCI_EXACT_FEATURE_DICTIONARY.csv`, phenotype-blind JSON |
| Literature pair eligibility and architecture encoding | `03_code/pairrank/pipeline.py::load_literature`, `internal_architecture`; `03_code/common/model.py::ArchitectureEncoder` | `validated_pairs_381.csv`, `eligible_non_tie_pairs_377.csv` |
| Scale-only ridge logistic architecture prior | `03_code/common/model.py::Standardizer`, `ArchitecturePrior`, `RidgeLogistic` | `overlap_exclusion_audit.csv` |
| Whole-system LOSO M0–M3 | `03_code/pairrank/pipeline.py::run_loso` | Internal16 pairs, PSCI matrix and system context |
| Mean-win-probability ranking and metrics | `03_code/pairrank/pipeline.py::rankings_and_metrics` | Frozen predictions and phenotype ranks |
| Computational check figures | `03_code/figures/render.py`; `04_scripts/05_reproduce_figures.py` | reproduced metrics and M2 probabilities |
| Acceptance comparison | `04_scripts/06_verify_reference_outputs.py`, `07_tests/` | model reference outputs |
