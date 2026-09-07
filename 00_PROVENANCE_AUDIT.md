# Frozen-repository provenance audit

Audit date: 2026-09-07. Repository locations below are logical labels, not server paths. Every listed commit exists in the shared Git history. Dirty status reflects only untracked material present at extraction; tracked source files at the stated commits were unchanged.

| repo | role | branch | frozen_commit | current_HEAD | commit_exists | clean_status | exact_source_files_used | notes |
|---|---|---|---|---|---|---|---|---|
| `historical_md` | original MD/descriptor executable lineage | `master` | `6fef6ba77965e504837e39cf2f5ef0205d86f567` | same | YES | TRACKED_CLEAN; unrelated untracked later material | `scripts/spatial/run_phase3b_spatial.py`; `scripts/spatial/analyze_spatial_metrics.py` | Defines executable D and L lineage and spatial extraction. |
| `truefit` | phenotype-blind PSCI resolution/Internal16 matrix | `psci-internal16-truefit-v1` | `c3fb6030860fa7d19bc0bb3782aa371487f24a93` | same | YES | CLEAN | `metadata/PSCI_PHENOTYPE_BLIND_FEATURE_FREEZE.json`; `scripts/psci_internal16_reanalysis.py`; `reports/PSCI_*_FEATURE_SPEC_FINAL.md`; `results/PSCI_FINAL/INTERNAL16_PSCI_MASTER_FINAL.csv` | Primary PSCI definition and value source. |
| `pairwise_v2` | frozen M0-M3 LOSO and release source | `system-conditioned-psci-v2` | `b895eb013f87f87da3f1e9a3e175bca77f2b30dd` | same | YES | TRACKED_CLEAN; untracked handoff material | `scripts/pairwise_v2/model.py`; `scripts/pairwise_v2/run_pairwise_v2.py`; `results/SYSTEM_CONDITIONED_PSCI_PAIRWISE_V2/*`; `results/NC_MANUSCRIPT_FREEZE_PATCH_V1/*` | `b895eb0` adds provenance/robustness audits without changing frozen M2. |
| `literature_v3` | Literature V3.2 final curation | `literature-v3-expansion` | `0e818da437056b0488a4033eac696beb91ada0ea` | same | YES | CLEAN | `data/literature_v3/construct_registry_v3_2_final.csv`; `phenotype_registry_raw_v3_2_final.csv`; `pairwise_labels_v3_2_final.csv` | Source for 91 constructs, 102 phenotype rows and 381 pairs. |
| `b_line` | final scientific-acceptance audit | `esm2-psci-rankreg-prep` | `f2f215ed4fa8b5ee86c36b2ec0c8fdbc8c965279` | same | YES | TRACKED_CLEAN; unrelated untracked handoff material | final acceptance reports referenced by the manuscript lineage | Audit anchor; not used to refit this release. |

The manuscript bundle was generated at Pairwise commit `b9ad43dfd7930824ecd0e2d452804d41e3bc4b23`; the release anchor is its descendant `b895eb0`. History is intentionally not rewritten or collapsed. Per-file mappings and source hashes are in `09_provenance/SOURCE_FILE_PROVENANCE.csv`.
