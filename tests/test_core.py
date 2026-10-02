import numpy as np
import pytest

from vec_metriclens.core import (
    contribution_rows,
    gene_attribution,
    metric_skill,
    spatial_compare,
)


def test_floor_and_ceiling_skill():
    assert metric_skill(0.0, 0.0, 0.8464, "higher") == pytest.approx(0.5)
    assert metric_skill(
        0.8464,
        0.0,
        0.8464,
        "higher",
    ) == pytest.approx(1.0)
    assert metric_skill(
        0.08359,
        0.08359,
        0.00406,
        "lower",
    ) == pytest.approx(0.5)


def test_t1_all_floor_metrics_sum_to_50():
    metrics = {
        "de_score": 0,
        "de_direction": 0,
        "mmd_u": 0.08359,
        "variogram": 0.005219,
    }
    rows, score = contribution_rows(metrics, "T1", "T1:val")
    assert score == pytest.approx(50.0)
    assert sum(row["points"] for row in rows) == pytest.approx(50.0)


def test_gene_attribution_identifies_largest_error():
    rows, summary = gene_attribution(
        np.array([3.0, 1.0]),
        np.array([2.0, 1.0]),
        np.array([1.0, 1.0]),
        ["a", "b"],
    )
    assert rows[0]["gene"] == "a"
    assert rows[0]["abs_change_error"] == pytest.approx(1.0)
    assert summary["change_mae"] == pytest.approx(0.5)


def test_spatial_scale_ratio_detects_doubling():
    rng = np.random.default_rng(2)
    truth = rng.normal(size=(300, 3))
    result = spatial_compare(truth * 2, truth)
    assert result["rms_scale_log_ratio"] == pytest.approx(
        np.log(2),
        abs=1e-10,
    )
    assert result["median_pairwise_distance_ratio"] == pytest.approx(
        2,
        rel=0.1,
    )
