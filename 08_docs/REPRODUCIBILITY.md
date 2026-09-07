# Reproducibility

Quick reproduction starts from frozen, trajectory-derived PSCI values and the frozen Literature V3.2 tables. It re-executes the specified prior, M0–M3 whole-system LOSO models, ranking aggregation, metrics, and compact computational figures. It does not tune or select a model.

The quick workflow is deterministic. Ridge logistic models have no intercept, use scale-only population SD and fixed L2 penalty 1.0. No random sampling occurs. The verification gate compares 96 M0–M3 pair probabilities, including all 24 M2 values, with absolute tolerance 1e-12.

Full trajectory reproduction is `NOT_AVAILABLE` in this primary release because no reduced trajectory set has been demonstrated to reproduce every frozen descriptor exactly. A separate deposit must provide all 16 standardized protein-only trajectories, analysis reference/topology files, component and active-site mapping, and file hashes. Until then, the full wrapper exits nonzero rather than substituting approximate inputs.
