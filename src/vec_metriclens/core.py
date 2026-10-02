from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

from .specs import ANCHORS, TASK_METRICS


def number(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def metric_skill(
    value: Any,
    floor: float,
    ceiling: float,
    direction: str,
) -> float:
    x = number(value)
    if x is None:
        return 0.0
    if direction == "zero":
        x, floor, ceiling = abs(x), abs(floor), abs(ceiling)
        direction = "lower"
    if direction == "higher":
        distance = max(0.0, ceiling - x)
        floor_distance = ceiling - floor
    elif direction == "lower":
        distance = max(0.0, x - ceiling)
        floor_distance = floor - ceiling
    else:
        raise ValueError(f"unknown direction: {direction}")
    if floor_distance <= 0:
        raise ValueError("invalid floor/ceiling anchors")
    return min(floor_distance / (floor_distance + distance), 1.0)


def contribution_rows(
    metrics: dict,
    task: str,
    board: str | None = None,
) -> tuple[list[dict], float | None]:
    task = task.upper()
    if task not in TASK_METRICS:
        raise ValueError("task must be T1/T2/T3")
    if board:
        if board not in ANCHORS:
            raise ValueError(f"unknown published validation board: {board}")
        if ANCHORS[board]["task"] != task:
            raise ValueError(f"{board} does not belong to {task}")

    rows: list[dict] = []
    total = 0.0
    for metric, (direction, weight, group) in TASK_METRICS[task].items():
        raw = number(metrics.get(metric))
        row = {
            "metric": metric,
            "group": group,
            "direction": direction,
            "weight": weight,
            "raw": raw,
        }
        if board:
            floor, ceiling = ANCHORS[board]["metrics"][metric]
            skill = metric_skill(raw, floor, ceiling, direction)
            points = 100 * weight * skill
            row.update(
                floor=floor,
                ceiling=ceiling,
                skill=skill,
                points=points,
                points_lost_vs_ceiling=100 * weight * (1 - skill),
            )
            total += points
        rows.append(row)
    return rows, total if board else None


def matrix_mean(matrix, chunk_rows: int = 256) -> np.ndarray:
    n_rows, n_cols = map(int, matrix.shape)
    if n_rows == 0:
        raise ValueError("matrix has zero cells")
    total = np.zeros(n_cols, dtype=np.float64)
    for start in range(0, n_rows, chunk_rows):
        part = matrix[start : start + chunk_rows]
        if sparse.issparse(part):
            total += np.asarray(part.sum(axis=0)).ravel()
        else:
            total += np.asarray(part, dtype=np.float64).sum(axis=0)
    return total / n_rows


def gene_attribution(
    pred_mean: np.ndarray,
    true_mean: np.ndarray,
    ref_mean: np.ndarray,
    genes: list[str],
) -> tuple[list[dict], dict]:
    if not (
        len(pred_mean) == len(true_mean) == len(ref_mean) == len(genes)
    ):
        raise ValueError("gene vectors must have the same length")

    true_change = true_mean - ref_mean
    pred_change = pred_mean - ref_mean
    error = pred_change - true_change
    rows = []
    for gene, tc, pc, err in zip(
        genes,
        true_change,
        pred_change,
        error,
        strict=True,
    ):
        same_direction = (
            bool(np.sign(tc) == np.sign(pc))
            if tc != 0
            else abs(pc) < 1e-12
        )
        rows.append(
            {
                "gene": gene,
                "true_change": float(tc),
                "pred_change": float(pc),
                "change_error": float(err),
                "abs_change_error": float(abs(err)),
                "abs_true_change": float(abs(tc)),
                "same_direction": same_direction,
            }
        )
    rows.sort(key=lambda x: (-x["abs_change_error"], x["gene"]))

    top = np.argsort(-np.abs(true_change))[: min(100, len(true_change))]
    direction_agreement = (
        float(
            np.mean(
                np.sign(true_change[top]) == np.sign(pred_change[top])
            )
        )
        if len(top)
        else float("nan")
    )
    if np.std(true_change) > 0 and np.std(pred_change) > 0:
        corr = float(np.corrcoef(true_change, pred_change)[0, 1])
    else:
        corr = float("nan")

    summary = {
        "change_mae": float(np.mean(np.abs(error))),
        "change_rmse": float(np.sqrt(np.mean(error**2))),
        "change_pearson": corr,
        "top100_true_shift_direction_agreement": direction_agreement,
    }
    return rows, summary


def spatial_summary(coords, seed: int = 0) -> dict:
    c = np.asarray(coords, dtype=np.float64)
    if c.ndim != 2 or c.shape[1] < 3 or len(c) < 2:
        raise ValueError("coordinates must be cells x >=3")
    c = c[:, :3]
    if not np.isfinite(c).all():
        raise ValueError("coordinates contain non-finite values")

    centered = c - c.mean(axis=0, keepdims=True)
    rms_radius = float(
        np.sqrt(np.mean(np.sum(centered**2, axis=1)))
    )
    covariance = np.cov(centered, rowvar=False)
    eigenvalues = np.sort(np.linalg.eigvalsh(covariance))[::-1]

    rng = np.random.default_rng(seed)
    pair_n = min(5000, max(100, len(c) * 5))
    a = rng.integers(0, len(c), size=pair_n)
    b = rng.integers(0, len(c), size=pair_n)
    keep = a != b
    distances = np.linalg.norm(c[a[keep]] - c[b[keep]], axis=1)
    quantiles = (
        np.quantile(distances, [0.1, 0.5, 0.9])
        if len(distances)
        else np.zeros(3)
    )
    return {
        "cells": len(c),
        "rms_radius": rms_radius,
        "covariance_eigenvalues": [float(x) for x in eigenvalues],
        "pairwise_distance_q10": float(quantiles[0]),
        "pairwise_distance_q50": float(quantiles[1]),
        "pairwise_distance_q90": float(quantiles[2]),
    }


def spatial_compare(pred, truth, seed: int = 0) -> dict:
    p = spatial_summary(pred, seed)
    t = spatial_summary(truth, seed)
    eps = 1e-12
    return {
        "prediction": p,
        "truth": t,
        "rms_scale_log_ratio": float(
            math.log(
                (p["rms_radius"] + eps) / (t["rms_radius"] + eps)
            )
        ),
        "median_pairwise_distance_ratio": float(
            (p["pairwise_distance_q50"] + eps)
            / (t["pairwise_distance_q50"] + eps)
        ),
    }


def write_rows(rows: list[dict], path: Path) -> None:
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if fields:
            writer.writeheader()
            writer.writerows(rows)
