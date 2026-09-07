# ACS manuscript availability text

## Code Availability — replace bracketed fields after public deposit

Custom code used to assemble PSCI descriptors from deposited descriptor inputs, fit the literature-derived architecture prior, perform whole-enzyme-system leave-one-system-out PSCI-PairRank evaluation, reproduce the reported computational metrics, and generate computational check figures is available at **[GITHUB_URL]**, release **v1.0.1-acs-reproducibility**, commit/release archive **[SOFTWARE_DOI]**. The archived release includes the pinned environment, exact run instructions, frozen model configuration, deterministic seed policy, tests, reference outputs, and SHA-256 manifests. In the pinned environment, the frozen probability outputs reproduce at the 1e-12 verification gate; compatible modern numerical stacks retain identical pair directions and manuscript metrics within the documented floating-point tolerance.

## Data Availability — compact release

The archived reproducibility release contains the Literature V3.2 architecture tables used for the prior, the Internal16 construct/rank freeze, experimental endpoint summaries, MD-QC summaries, the exact PSCI feature dictionary and 16-topology PSCI matrix, all 24 frozen M2 out-of-fold pair predictions, M0–M3 comparison tables, fold coefficients/scalers, current Figure 1–6 computational source tables, provenance mappings and cryptographic manifests. Individual experimental replicate measurements and SDs are included only when recoverable from authoritative primary records and are not reconstructed from summary means. The compact release begins the quick reproduction from frozen trajectory-derived PSCI descriptors. The 16 standardized MD trajectories required for trajectory-to-descriptor recreation are **[TRAJECTORY_DOI / TO BE DEPOSITED SEPARATELY]**.

## Before final publication

Replace `[GITHUB_URL]`, `[SOFTWARE_DOI]`, `[TRAJECTORY_DOI]` and the related article DOI in `CITATION.cff`/`zenodo.json`. If the trajectory archive is not deposited at initial submission, do not describe the compact package as full trajectory-to-feature reproducibility.
