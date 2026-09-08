# Final public-release audit

## Release scope

- `VERSION = 1.0.2`
- `RELEASE_TYPE = DOCUMENTATION_ONLY_MAINTENANCE`
- `SOFTWARE_CREATORS = Cai, Xue; Xie, Chengqian`
- `SCIENTIFIC_NUMERIC_FILES_CHANGED = 0`
- `MODEL_FILES_CHANGED = 0`
- `CODE_LOGIC_FILES_CHANGED = 0`
- `CREATOR_ATTRIBUTION_CONSISTENCY = EXACT`
- `CURRENT_SOFTWARE_ATTRIBUTION_CONFLICTS = 0`
- `FULL_TRAJECTORY_REPRODUCTION = NOT_AVAILABLE`
- `ARTICLE_DOI_STATUS = NOT_YET_AVAILABLE`
- `TRAJECTORY_DOI_STATUS = NOT_YET_AVAILABLE`

## Frozen scientific anchors

- Internal16: four systems, 16 topologies and 24 pairs.
- Literature V3.2: 19 studies, 20 systems, 91 constructs, 102 phenotype records, 381 validated pairs, 377 eligible non-ties and four ties.
- M2: 19/24 pair accuracy; Spearman 0.6548461876; Kendall tau-b 0.5217491947; Top1 3/4; true-best in predicted Top2 3/4; exact Top2 set 2/4.

Pinned verification uses Python 3.12.3, NumPy 1.26.4 and SciPy 1.11.4 with probability tolerance `1e-12`. Compatible modern environments use `1e-7` and must retain identical pair directions and manuscript metrics.
