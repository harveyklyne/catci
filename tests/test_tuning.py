"""The CV tuner (catci.tuning) and the d axis of the experiment configs."""

import numpy as np
import pytest

xgb = pytest.importorskip("xgboost")

import dgp
from catci.tuning import _stack, kfold_indices, tune_xgboost
from config import Config, d_grid, power_config, size_config, tuning_path

TINY_GRID = dict(eta=[0.3], max_depth=[1, 2], gamma=[0.0])


def test_kfold_indices_partition():
    folds = kfold_indices(103, 5, np.random.default_rng(0))
    tests = np.concatenate([te for _, te in folds])
    assert sorted(tests.tolist()) == list(range(103))
    for tr, te in folds:
        assert len(np.intersect1d(tr, te)) == 0 and len(tr) + len(te) == 103


def test_stacked_folds_stay_inside_their_replicate():
    # each booster must only ever see one replicate: that is what keeps the
    # training size at (K-1)/K * n rather than growing with the number of replicates
    rng = np.random.default_rng(1)
    sizes = [40, 55, 30]
    data = [(rng.normal(size=(m, 5)), rng.integers(1, 4, size=m)) for m in sizes]
    z, labels, folds = _stack(data, 5, rng)
    assert z.shape == (sum(sizes), 5) and len(folds) == 5 * len(sizes)
    bounds = np.cumsum([0] + sizes)
    for r in range(len(sizes)):
        for tr, te in folds[5 * r: 5 * (r + 1)]:
            both = np.concatenate([tr, te])
            assert both.min() >= bounds[r] and both.max() < bounds[r + 1]
            assert len(both) == sizes[r]


def test_tune_xgboost_picks_from_grid_and_records_table():
    rng = np.random.default_rng(2)
    datasets = []
    for _ in range(2):
        sim = dgp.simulate_marginal(400, 4, "lin", rng)  # lin: Z moves the pmf a lot
        datasets.append((sim["z"], sim["x"]))
    res = tune_xgboost(datasets, 4, grid=TINY_GRID, max_rounds=60, early_stopping_rounds=10, rng=rng)
    assert len(res.table) == 2
    assert res.params["eta"] == 0.3 and res.params["max_depth"] in (1, 2)
    assert 1 < res.params["nrounds"] <= 60
    assert res.cv_logloss == min(row["logloss"] for row in res.table)
    # beats the Z-blind predictor (the marginal class frequencies)
    x = np.concatenate([x for _, x in datasets])
    freq = np.bincount(x, minlength=5)[1:] / len(x)
    assert res.cv_logloss < -np.sum(freq * np.log(freq))
    js = res.to_json()
    assert set(js["xgb"]) == {"eta", "max.depth", "gamma", "nrounds"}


def test_tune_xgboost_tolerates_classes_missing_from_a_fold():
    # at d in the hundreds many classes have fewer members than folds
    rng = np.random.default_rng(3)
    z = rng.normal(size=(60, 5))
    labels = rng.integers(1, 6, size=60)
    labels[:2] = 9  # class 9 appears twice, classes 6-8 and 10 never
    res = tune_xgboost((z, labels), 10, grid=dict(eta=[0.3], max_depth=[1], gamma=[0.0]),
                       max_rounds=20, early_stopping_rounds=None, rng=rng)
    assert np.isfinite(res.cv_logloss)


def test_tune_xgboost_rejects_out_of_range_labels():
    with pytest.raises(ValueError):
        tune_xgboost((np.zeros((10, 1)), np.arange(10)), 10, grid=TINY_GRID, max_rounds=5)


# --------------------------------------------------------------------------- #
# the d axis of the configs
# --------------------------------------------------------------------------- #
def test_default_configs_keep_their_names_and_tuning():
    cfg = power_config("sin", "lin", "binary_tree")
    assert (cfg.n, cfg.dx, cfg.dy) == (1000, 8, 8)
    assert cfg.name == "power_sin_lin_binary_tree"
    assert size_config("lin", "lin").name == "size_lin_lin"
    # the frozen R tuning is still what the default reads
    assert cfg.xgb_params("sin", 8) == {"eta": 0.01, "max.depth": 1, "gamma": 1.5, "nrounds": 861}


def test_d_grid_varies_dx_only():
    cfgs = d_grid("lin", "lin", "step", [8, 32, 128], dy=4, n=2000, reps=10)
    assert [(c.n, c.dx, c.dy, c.reps) for c in cfgs] == [(2000, 8, 4, 10), (2000, 32, 4, 10), (2000, 128, 4, 10)]
    assert len({c.name for c in cfgs}) == 3
    assert cfgs[1].name == "power_lin_lin_step_n2000_dx32_dy4"


def test_missing_tuning_names_the_command():
    cfg = power_config("lin", "lin", "step", n=12345, d=8)
    assert not tuning_path(12345, 8, "lin").exists()
    with pytest.raises(FileNotFoundError, match="tune.py --n 12345 --d 8"):
        cfg.xgb_params("lin", 8)
