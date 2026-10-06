import subprocess
import sys

import pytest

ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")


def _write(path, expression, genes, labels=None):
    obs = pd.DataFrame(
        index=[f"cell_{i}" for i in range(len(expression))]
    )
    if labels is not None:
        obs["celltype"] = labels
    ad.AnnData(
        X=np.asarray(expression, dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=genes),
    ).write_h5ad(path)


def test_cli_runs_real_scorer_and_writes_gene_lens(tmp_path):
    rng = np.random.default_rng(11)
    cells, n_genes = 100, 32
    genes = [f"g{i}" for i in range(n_genes)]
    labels = np.array(["A"] * 50 + ["B"] * 50)

    reference = np.log1p(
        rng.poisson(3, size=(cells, n_genes))
    ).astype(np.float32)
    target = reference.copy()
    target[:, :8] += 0.45
    prediction = np.maximum(
        target + rng.normal(0, 0.1, target.shape),
        0,
    )

    ref = tmp_path / "ref.h5ad"
    truth = tmp_path / "target.h5ad"
    pred = tmp_path / "pred.h5ad"
    _write(ref, reference, genes, labels)
    _write(truth, target, genes, labels)
    _write(pred, prediction, genes)

    out = tmp_path / "out"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "vec_metriclens.cli",
            str(pred),
            "--task",
            "T1",
            "--target",
            str(truth),
            "--reference",
            str(ref),
            "--top-genes",
            "5",
            "--out",
            str(out),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )

    assert run.returncode == 0, run.stdout + run.stderr
    assert (out / "genes.csv").exists()
    report = (out / "report.md").read_text()
    assert "Gene-change attribution" in report
    assert "Raw ranking metrics" in report



def test_cli_rejects_t2_board_setting_mismatch(tmp_path):
    dummy = tmp_path / "dummy.h5ad"
    ref = tmp_path / "ref.h5ad"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "vec_metriclens.cli",
            str(dummy),
            "--task",
            "T2",
            "--setting",
            "embryo",
            "--target",
            str(dummy),
            "--reference",
            str(ref),
            "--board",
            "T2:heart:val_extrap",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode != 0
    assert "belongs to T2 setting" in (run.stdout + run.stderr)
