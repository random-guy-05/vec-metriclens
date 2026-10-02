# Community Contribution submission text

## Title
VEC MetricLens — explain local scores with point attribution, gene errors and spatial diagnostics

## Description
VEC MetricLens runs the real public `veckit==0.1.2` scorer and turns its metric panel into an actionable debugging report. For any pseudo-target it labels each ranking metric by its official biological question, optimization direction and weight, then computes per-gene predicted-vs-observed change against the same reference/WT and ranks the largest gene-level errors. For T2/T3 it additionally reports translation-invariant spatial scale, covariance-eigenvalue and pairwise-distance diagnostics. When—and only when—the supplied target is a published validation board, `--board` reproduces the public hyperbolic skill transform and official task weights to show each metric's point contribution and points lost versus the published ceiling. It deliberately refuses to invent an aggregate for unrelated pseudo-targets. The tool helps entrants move from “the score changed” to “which biological/scoring component is limiting this model, and where should I debug next?”
