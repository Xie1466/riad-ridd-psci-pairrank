# Changelog

## v1.0.2 — 2026-09-08 documentation-only maintenance release

- Unified software creator metadata to the authoritative Zenodo order: Cai, Xue; Xie, Chengqian.
- Updated version, citation, persistent-identifier, reproducibility-tolerance and PSCI-definition-precedence documentation.
- No scientific data, model parameters, PSCI values, pair labels, predictions, reported metrics or executable model logic changed.

## v1.0.1 — 2026-09-08 DOI-frozen documentation-only cleanup

- Updated publication-state wording to the issued Zenodo DOI without changing the GitHub tag, commit, Zenodo files, scientific data, model code, or results.
- Synchronized `.zenodo.json` to the published Zenodo description, keywords, and DOI creator list.
- Clarified DOI-facing versus repository/software attribution, numerical tolerance modes, final PSCI definition priority, and missingness policy. **SUPERSEDED BY V1.0.2 SOFTWARE CREATOR POLICY.**

## v1.0.1 — 2026-09-07 public-release audit patch

- No model, feature, penalty, split, phenotype, PSCI value or frozen prediction was changed.
- Synchronized `05_data/figures/` to the current ACS manuscript Figure 1–6.
- Removed historical server absolute paths from the stored reference PSCI output.
- Added author-verified experimental reporting annotations (`n = 3`, common `(G4S)3` experimental linker) without fabricating raw replicates or SDs.
- Added explicit data-license boundaries and final author metadata for software citation.
- Added exact-runtime detection and a portability gate: strict 1e-12 under the pinned runtime, 1e-7 probability tolerance under alternate modern numerical stacks while pair directions/metrics must remain identical.
- Removed a stale reference to a non-bundled clean-room report and added `PORTABILITY_AUDIT.md`.

## 1.0.0 — 2026-09-07

- Initial ACS reproducibility freeze.
- Includes the frozen PSCI matrix, Literature V3.2 architecture-prior inputs, Internal16 M0–M3 out-of-fold references, source data, public wrappers, tests, and provenance.
- Excludes production trajectories pending a separate DOI-backed trajectory deposit.
