# Community Contribution submission text

## Title
VEC MetricLens — gene-level and spatial debugging for local VEC scores

## Description
VEC MetricLens runs the real public `veckit==0.1.2` scorer and turns weak local metrics into concrete places to debug in the prediction itself. It computes every gene's predicted and observed pseudobulk change against the same reference/WT, ranks the largest gene-level errors, summarizes change correlation and direction agreement, and for T2/T3 reports translation-invariant spatial scale, covariance-eigenvalue and pairwise-distance diagnostics. It also labels each official metric by its biological question, direction and weight. When—and only when—the supplied target is a published validation board, `--board` reproduces the public skill transform and point accounting; unrelated pseudo-targets get no invented aggregate. This is distinct from comparing two scorer JSON files: its main output is attribution inside one prediction, helping entrants identify which genes or spatial properties are responsible for a weak score.
