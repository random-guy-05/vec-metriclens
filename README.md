# VEC MetricLens

**Turn a local VEC score into a debugging report: which scoring questions cost points, which genes are wrong, and—when spatial—which geometry differs.**

MetricLens runs the real public `veckit==0.1.2` scorer and adds three layers of interpretation that the raw JSON does not provide:

1. **ranking-metric lens** — raw metric, direction, official weight/question, and optionally exact published validation-board skill/point contribution;
2. **gene-change attribution** — predicted vs true pseudobulk change against the same reference/WT, ranked by gene-level change error;
3. **spatial geometry lens** for T2/T3 — translation-invariant RMS radius, covariance eigenvalues and sampled pairwise-distance scale.

These diagnostics do not create new competition metrics. They explain discrepancies in the local files you supplied.

## Usage

```bash
pip install -e .

# Any local pseudo-target: raw metrics + attribution, no fake aggregate
vec-metriclens prediction.h5ad   --task T1   --target pseudo_target.h5ad   --reference preceding_stage.h5ad   --out lens

# Published validation-board accounting (only if the supplied target really is that board)
vec-metriclens prediction.h5ad   --task T2 --setting heart   --target validation_target.h5ad   --reference preceding_stage.h5ad   --board T2:heart:val_extrap   --out lens
```

For Task 3 use `--wt matched_wt.h5ad`.

Outputs:

- `analysis.json` — scorer metrics + structured diagnostics;
- `genes.csv` — every gene's true change, predicted change, error and direction agreement;
- `report.md` — ranked bottlenecks and top gene/spatial discrepancies.

## Published-board contribution accounting

When `--board` is supplied, MetricLens reproduces the current published hyperbolic skill transform and official task weights. The report sorts metrics by **points lost versus the published ceiling**, which turns “my score is 61” into an actionable bottleneck list.

Do **not** apply a validation board's anchors to an unrelated pseudo-target. Without `--board`, MetricLens deliberately refuses to invent a combined score and reports only raw local metrics.

## Gene lens

For each gene:

- true change = mean(target) − mean(reference/WT);
- predicted change = mean(prediction) − mean(reference/WT);
- change error = predicted − true change;
- sign agreement tells whether the response moved in the same direction.

The report also gives global change MAE/RMSE, predicted-vs-true change correlation, and direction agreement among the 100 strongest observed shifts. Large H5AD matrices are mean-reduced in chunks instead of densifying the full transcriptome at once.

## Spatial lens

For T2/T3, geometry summaries use only the first three `spatial_3D` columns and center each cloud, so arbitrary translations do not contaminate the diagnostics. The lens reports overall scale and anisotropy-like eigenvalue structure; the official scorer remains the authority for the actual spatial metrics.

See `docs/SOURCES.md` for the official metric/anchor snapshot.
