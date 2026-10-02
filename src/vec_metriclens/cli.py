from __future__ import annotations

import argparse
import json
from contextlib import suppress
from pathlib import Path

from .core import (
    contribution_rows,
    gene_attribution,
    matrix_mean,
    spatial_compare,
    write_rows,
)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Explain a local VEC score beyond one aggregate number."
    )
    p.add_argument("prediction", type=Path)
    p.add_argument("--task", required=True, choices=["T1", "T2", "T3"])
    p.add_argument("--setting", default="heart", choices=["heart", "embryo"])
    p.add_argument("--target", required=True, type=Path)
    p.add_argument("--reference", type=Path)
    p.add_argument("--wt", type=Path)
    p.add_argument(
        "--board",
        help=(
            "published validation board; only use when the supplied target "
            "is that board"
        ),
    )
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--top-genes", type=int, default=20)
    p.add_argument("--out", type=Path, default=Path("metriclens_out"))
    return p


def _open(path: Path):
    import anndata as ad

    return ad.read_h5ad(path, backed="r")


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.task in {"T1", "T2"} and args.reference is None:
        raise SystemExit("--reference is required for T1/T2")
    if args.task == "T3" and args.wt is None:
        raise SystemExit("--wt is required for T3")
    if args.top_genes < 1:
        raise SystemExit("--top-genes must be >= 1")

    from veckit import score

    kwargs = {
        "task": args.task,
        "input": args.prediction,
        "target": args.target,
        "seed": args.seed,
    }
    ref_path = args.reference if args.task in {"T1", "T2"} else args.wt
    if args.task in {"T1", "T2"}:
        kwargs["reference"] = args.reference
    else:
        kwargs["wt"] = args.wt
    if args.task == "T2":
        kwargs["setting"] = args.setting

    result = score(**kwargs)
    metrics = result["metrics"]
    contribution, official_score = contribution_rows(
        metrics,
        args.task,
        args.board,
    )

    pred = _open(args.prediction)
    truth = _open(args.target)
    ref = _open(ref_path)
    try:
        genes = [str(x) for x in pred.var_names]
        truth_genes = [str(x) for x in truth.var_names]
        ref_genes = [str(x) for x in ref.var_names]
        if genes != truth_genes or genes != ref_genes:
            raise ValueError(
                "prediction/target/reference gene names and order "
                "must match for attribution"
            )

        pred_mean = matrix_mean(pred.X)
        true_mean = matrix_mean(truth.X)
        ref_mean = matrix_mean(ref.X)
        gene_rows, gene_summary = gene_attribution(
            pred_mean,
            true_mean,
            ref_mean,
            genes,
        )

        spatial = None
        if args.task in {"T2", "T3"}:
            if (
                "spatial_3D" not in pred.obsm
                or "spatial_3D" not in truth.obsm
            ):
                raise ValueError(
                    "T2/T3 attribution requires obsm['spatial_3D'] "
                    "in prediction and target"
                )
            spatial = spatial_compare(
                pred.obsm["spatial_3D"],
                truth.obsm["spatial_3D"],
                args.seed,
            )
    finally:
        for data in (pred, truth, ref):
            with suppress(AttributeError, OSError, ValueError):
                data.file.close()

    args.out.mkdir(parents=True, exist_ok=True)
    write_rows(gene_rows, args.out / "genes.csv")

    payload = {
        "meta": result.get("meta", {}),
        "metrics": metrics,
        "metric_contributions": contribution,
        "published_validation_score": official_score,
        "gene_summary": gene_summary,
        "spatial_summary": spatial,
    }
    (args.out / "analysis.json").write_text(
        json.dumps(payload, indent=2, default=float) + "\n",
        encoding="utf-8",
    )

    ranked_contribution = sorted(
        contribution,
        key=lambda row: row.get("points_lost_vs_ceiling", -1),
        reverse=True,
    )
    lines = ["# VEC MetricLens", ""]
    if official_score is not None:
        lines.extend(
            [
                (
                    "Published validation-scale score: "
                    f"**{official_score:.3f}**"
                ),
                "",
                (
                    "> This transform is board-specific. Do not apply "
                    "validation anchors to an unrelated pseudo-target."
                ),
                "",
                "## Ranking-metric contribution",
                "",
                (
                    "| metric | question | raw | skill | points | "
                    "points lost vs ceiling |"
                ),
                "|---|---|---:|---:|---:|---:|",
            ]
        )
        for row in ranked_contribution:
            raw = "NA" if row["raw"] is None else f"{row['raw']:.6g}"
            lines.append(
                f"| {row['metric']} | {row['group']} | {raw} | "
                f"{row['skill']:.3f} | {row['points']:.2f} | "
                f"{row['points_lost_vs_ceiling']:.2f} |"
            )
    else:
        lines.extend(
            [
                "## Raw ranking metrics",
                "",
                "| metric | question | direction | raw | weight |",
                "|---|---|---|---:|---:|",
            ]
        )
        for row in contribution:
            raw = "NA" if row["raw"] is None else f"{row['raw']:.6g}"
            lines.append(
                f"| {row['metric']} | {row['group']} | "
                f"{row['direction']} | {raw} | {row['weight']:.3f} |"
            )

    lines.extend(
        [
            "",
            "## Gene-change attribution",
            "",
            (
                "- mean absolute change error: "
                f"**{gene_summary['change_mae']:.6g}**"
            ),
            (
                "- RMSE of gene-level change: "
                f"**{gene_summary['change_rmse']:.6g}**"
            ),
            (
                "- Pearson correlation of predicted vs true gene change: "
                f"**{gene_summary['change_pearson']:.4f}**"
            ),
            (
                "- direction agreement among 100 strongest true shifts: "
                f"**{gene_summary['top100_true_shift_direction_agreement']:.1%}**"
            ),
            "",
            f"### Top {args.top_genes} genes by absolute change error",
            "",
            (
                "| gene | true change | predicted change | "
                "absolute error | direction correct |"
            ),
            "|---|---:|---:|---:|---|",
        ]
    )
    for row in gene_rows[: args.top_genes]:
        lines.append(
            f"| {row['gene']} | {row['true_change']:.5g} | "
            f"{row['pred_change']:.5g} | "
            f"{row['abs_change_error']:.5g} | "
            f"{'yes' if row['same_direction'] else 'no'} |"
        )

    if spatial:
        pred_spatial = spatial["prediction"]
        truth_spatial = spatial["truth"]
        lines.extend(
            [
                "",
                "## Spatial geometry lens",
                "",
                (
                    "- prediction RMS radius: "
                    f"**{pred_spatial['rms_radius']:.5g}**; target: "
                    f"**{truth_spatial['rms_radius']:.5g}**"
                ),
                (
                    "- RMS scale log-ratio: "
                    f"**{spatial['rms_scale_log_ratio']:.5g}** "
                    "(0 is scale-matched)"
                ),
                (
                    "- median pairwise-distance ratio: "
                    f"**{spatial['median_pairwise_distance_ratio']:.5g}**"
                ),
                (
                    "- prediction covariance eigenvalues: "
                    f"**{pred_spatial['covariance_eigenvalues']}**"
                ),
                (
                    "- target covariance eigenvalues: "
                    f"**{truth_spatial['covariance_eigenvalues']}**"
                ),
            ]
        )

    lines.extend(
        [
            "",
            (
                "> MetricLens is diagnostic. Gene and geometry summaries "
                "explain local discrepancies; they are not additional "
                "official metrics and do not expose hidden targets."
            ),
            "",
        ]
    )
    (args.out / "report.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
