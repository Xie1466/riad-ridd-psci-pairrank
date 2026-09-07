# Changelog

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