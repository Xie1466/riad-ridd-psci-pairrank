# RIAD/RIDD PSCI-PairRank ACS reproducibility release v1.0.2

This is the compact code-and-data release for the computational component of the RIAD/RIDD protein-scaffold study **“Productive spatial coupling shapes topology-dependent performance of modular enzyme scaffolds.”** It contains Literature V3.2 (19 studies, 20 systems, 91 constructs, 102 phenotype records and 381 validated pairs), the Internal16 four-system/16-topology freeze, the structured PSCI representation, the frozen literature architecture prior and the final M2 PSCI-PairRank evaluation.

## Current public release

- GitHub repository: https://github.com/Xie1466/riad-ridd-psci-pairrank
- Release: `v1.0.2-acs-reproducibility`
- Release commit: `70fb08f246f88302c6797b5ebccc9774a409637b`
- GitHub release: https://github.com/Xie1466/riad-ridd-psci-pairrank/releases/tag/v1.0.2-acs-reproducibility
- Zenodo record: https://zenodo.org/records/22656378
- Zenodo version DOI: https://doi.org/10.5281/zenodo.22656378
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.22638943
- Software creators: Cai, Xue; Xie, Chengqian

The previous release, `v1.0.1-acs-reproducibility` (https://doi.org/10.5281/zenodo.22638944), is superseded documentation-wise and remains a scientifically equivalent frozen predecessor. Version 1.0.2 changes documentation, citation attribution, persistent-identifier metadata and release consistency only.

## What this release reproduces

The quick workflow starts from frozen trajectory-derived PSCI values and public literature-architecture tables and re-executes PSCI feature-table assembly and validation, the literature architecture prior, M0–M3 whole-enzyme-system LOSO evaluation, construct-level mean-win-probability ranking, reported computational metrics and compact computational check figures.

Expected M2 results are **19/24 pair directions (79.2%)**, Spearman **0.6548461876**, Kendall tau-b **0.5217491947**, Top1 **3/4**, true-best-in-predicted-Top2 **3/4**, and exact Top2 set **2/4**.

## Recommended exact reproduction

Use the pinned environment whenever possible:

```bash
conda env create -f 01_environment/environment.yml
conda activate riad-ridd-psci-reproduction
bash 04_scripts/run_quick_reproduction.sh
PYTHONPATH=07_tests python -m unittest discover -s 07_tests -p 'test_*.py'
```

A `venv` is acceptable only with Python 3.12.3. The pinned Python 3.12.3, NumPy 1.26.4 and SciPy 1.11.4 environment uses the strict `1e-12` probability gate. Compatible modern numerical stacks use `1e-7` while still requiring identical pair directions and manuscript metrics; see `08_docs/PORTABILITY_AUDIT.md`.

## Scientific definition authority

The final authoritative PSCI definitions are `05_data/internal16/PSCI_EXACT_FEATURE_DICTIONARY.csv` and `03_code/psci/descriptors.py`. The phenotype-blind JSON in `02_config/` is historical provenance and cannot override the final executable definition. In particular, `D` is the mass-weighted enzyme all-atom center-of-mass distance normalized by reference Cα radius of gyration, and `L` is linker end-to-end distance normalized by instantaneous summed Cα–Cα contour length. The predictive PSCI feature set is `[D, C, F1, F4, F5, L]`; F3 remains not uniquely recoverable, Q_topology is not fabricated, and A/O missingness remains evidence-qualified.

## Data scope

The compact release does not synthesize missing replicate values or uncertainty statistics. Production trajectories are not included; full trajectory-to-feature reproduction is `NOT_AVAILABLE`. The required trajectory/reference/topology/index inventory is documented under `05_data/md_reduced_or_manifest/`.

Software is MIT licensed. Project-authored derived data are CC BY 4.0; literature-derived metadata retain source-specific licensing terms described in `DATA_LICENSE.md`. `ARTICLE_DOI_STATUS = NOT_YET_AVAILABLE` and `TRAJECTORY_DOI_STATUS = NOT_YET_AVAILABLE`.
