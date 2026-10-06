# MetricLens ranking contract

Snapshot date: **2026-10-05**.

MetricLens does not redefine the VEC metrics. `veckit==0.1.2` computes them; MetricLens only labels, attributes, and (when a real published validation board is supplied) converts them through the public skill map.

## Task 1

| Metric | Direction | Global weight | Question |
|---|---|---:|---|
| `de_score` | higher | 0.25 | DE gene recovery |
| `de_direction` | higher | 0.25 | Change direction |
| `mmd_u` | lower | 0.30 | Cell-state distribution |
| `variogram` | lower | 0.20 | Gene-gene co-variation |

## Task 2

| Metric | Direction | Global weight | Question |
|---|---|---:|---|
| `de_score` | higher | 0.125 | Expression change |
| `de_direction` | higher | 0.125 | Expression change |
| `mmd_u` | lower | 0.15 | Cell-state distribution |
| `variogram` | lower | 0.10 | Cell-state distribution group |
| `d2_shape` | lower | 1/12 | Tissue shape and growth scale |
| `occupancy_dice` | higher | 1/12 | Tissue shape and growth scale |
| `scale_log_ratio` | zero | 1/12 | Tissue shape and growth scale |
| `neighborhood_mmd` | lower | 0.25 | Local spatial organisation |

## Task 3

| Metric | Direction | Global weight | Question |
|---|---|---:|---|
| `de_score` | higher | 0.30 | Response gene recovery |
| `de_direction` | higher | 0.25 | Response direction |
| `severity_slope` | zero | 0.25 | Response magnitude |
| `mmd_u` | lower | 0.12 | Cell-state distribution |
| `variogram` | lower | 0.08 | Cell-state distribution group |

## Published validation skill map

The official score first maps each raw metric to a bounded skill using the board-specific floor and ceiling, then applies the task weights. Floor performance maps to skill 0.5; ceiling maps to 1.0. Signed zero-ideal terms are converted to absolute magnitude before the lower-is-better transform.

MetricLens stores the published validation anchors in `src/vec_metriclens/specs.py`. Hidden test anchors are not known and are not guessed.

## MetricLens-only diagnostics

The following are explanatory summaries, **not official metrics**:

- gene-level pseudobulk change error;
- top-100 shift direction agreement;
- change-vector Pearson correlation;
- centered RMS spatial radius;
- covariance eigenvalues;
- sampled pairwise-distance ratio.

They are intended to help locate failure modes after the official metric panel has identified a weak question.
