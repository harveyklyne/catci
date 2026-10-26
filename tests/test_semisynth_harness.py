"""The semi-synthetic harness after the random baselines were added (RANDOM_SEARCH.md).

* the report's paired baseline is the full-depth merge, not the first ``merge@`` column;
* the sweep's existing methods are bit-for-bit what they were before ``merge@k`` /
  ``rsplit@k`` / ``rmerge@k`` were added (golden replicate, captured from the code at
  the parent commit);
* ``spawn(4)[:3] == spawn(3)``, the fact that golden replicate relies on.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

import adult_semisynth as ss
import report_semisynth as rs
import sweep_semisynth as sw


def test_report_baseline_is_full_depth_merge():
    cols = ["split@1", "merge@1", "merge@2", "merge@10", "merge@14", "rmerge@1", "max@0"]
    assert rs.full_merge(cols) == "merge@14"  # not merge@1 (first) nor merge@2 (alphabetical last)
    assert rs.full_merge(["merge@14"]) == "merge@14"


def test_route_differences_pair_within_route():
    rng = np.random.default_rng(0)
    wide = pd.DataFrame({c: rng.integers(0, 2, 50).astype(float)
                         for c in ["split@1", "rsplit@1", "merge@1", "rmerge@1", "merge@5"]})
    out = rs.route_differences(wide)
    assert list(out.index) == ["split@1 - rsplit@1", "merge@1 - rmerge@1"]  # merge@5 has no random twin
    d = wide["split@1"] - wide["rsplit@1"]
    assert np.isclose(out.loc["split@1 - rsplit@1", "diff"], d.mean())
    assert np.isclose(out.loc["split@1 - rsplit@1", "se"], d.std(ddof=1) / np.sqrt(50))


# Education x Income, planted, n = 300, lam = 0.8, depths 1 2, B in (50, 100),
# SeedSequence(123), replicate 0. Captured by running experiments/sweep_semisynth.py
# as of commit 6dca08e (before the random arms) through the same one_rep call.
GOLDEN = {
    ("split", 1, 50): 0.0196078431372549, ("split", 1, 100): 0.009900990099009901,
    ("split", 2, 50): 0.0196078431372549, ("split", 2, 100): 0.009900990099009901,
    ("split", 14, 50): 0.0392156862745098, ("split", 14, 100): 0.04950495049504951,
    ("merge", 14, 50): 0.0392156862745098, ("merge", 14, 100): 0.039603960396039604,
    ("max", 0, 50): 0.13725490196078433, ("max", 0, 100): 0.13861386138613863,
    ("euclid", 0, 50): 0.17647058823529413, ("euclid", 0, 100): 0.18811881188118812,
    ("mGCM", 0, 50): 0.3333333333333333, ("mGCM", 0, 100): 0.3564356435643564,
}


def test_sweep_existing_methods_unchanged_and_new_methods_present():
    args = argparse.Namespace(x="Education", y="Income", w=list(ss.DEFAULT_W), z=list(ss.DEFAULT_W),
                              direction="planted", n=300, pool=5.0, depths=[1, 2],
                              n_boots=[50, 100], time_reps=1, taxonomy=True)
    rows, timing = sw.one_rep((args, 0.8, 0, np.random.SeedSequence(123)))
    got = {(r["method"], r["k"], r["n_boot"]): r["p_value"] for r in rows}
    for key, p in GOLDEN.items():
        assert got[key] == p, key  # exact: same tie-break streams, same paths
    methods = {(r["method"], r["k"]) for r in rows}
    L = 14
    assert {("merge", 1), ("merge", 2), ("rsplit", 1), ("rsplit", L), ("rmerge", 2), ("rmerge", L)} <= methods
    assert not any(m == "mergecoarse" for m, _ in methods)
    assert ("merge", L) in methods and sum(m == "merge" and k == L for m, k in methods) == 1
    assert all(0.0 < r["p_value"] <= 1.0 for r in rows)
    # full searches carry their own depth; direct runs exist for the truncated merge and random arms
    assert {t["k"] for t in timing if t["method"] in ("split", "merge", "rsplit", "rmerge")} == {L}
    direct = {(t["method"], t["k"]) for t in timing if t["method"].endswith("_direct")}
    assert {("merge_direct", 1), ("rsplit_direct", 2), ("rmerge_direct", L)} <= direct


def test_spawn4_children_0_to_2_equal_spawn3():
    for seed in (np.random.SeedSequence(7), np.random.SeedSequence(7, spawn_key=(3,))):
        three = np.random.default_rng(seed).spawn(3)
        four = np.random.default_rng(
            np.random.SeedSequence(seed.entropy, spawn_key=seed.spawn_key)).spawn(4)
        for a, b in zip(three, four[:3]):
            assert (a.bit_generator.seed_seq.spawn_key == b.bit_generator.seed_seq.spawn_key
                    and a.integers(2**62, size=5).tolist() == b.integers(2**62, size=5).tolist())
        assert four[3].bit_generator.seed_seq.spawn_key[-1] == 3


# run_semisynth.one_rep, oracle learner, Education x Income planted, n = 300, B = 100,
# lam = 0.8, SeedSequence(5), truncations 2 4; captured at commit 6dca08e likewise.
GOLDEN_RUN = {
    "merge": 0.13861386138613863, "split": 0.18811881188118812, "split@2": 0.3069306930693069,
    "split@4": 0.25742574257425743, "max": 0.2079207920792079, "euclid": 0.06930693069306931,
    "mGCM": 0.19801980198019803, "chi_sq": 0.06597256219625658, "at_typed": 0.0941258289429501,
    "at_cat": 0.065972562196238, "strat_chi2": 0.11605229348948587,
}


def test_run_semisynth_existing_methods_unchanged_and_random_block_present():
    import run_semisynth as rr
    args = argparse.Namespace(x="Education", y="Income", w=list(ss.DEFAULT_W), z=list(ss.DEFAULT_W),
                              direction="planted", n=300, pool=5.0, learners=["oracle", "oracle"],
                              truncations=[2, 4], full_split=True, n_boot=100, taxonomy=True,
                              random=True)
    rows = rr.one_rep((args, 0.8, 0, np.random.SeedSequence(5)))
    oracle = [(r["method"], r["p_value"]) for r in rows if r["learner"] == "oracle"]
    half = len(oracle) // 2
    first = dict(oracle[:half], strat_chi2=rows[-1]["p_value"])
    for name, p in GOLDEN_RUN.items():
        assert first[name] == p, name
    # the same learner twice sees identical random paths and tie-breaks
    assert oracle[:half] == oracle[half:]
    assert {"merge@2", "merge@4", "rsplit@2", "rsplit@4", "rmerge@2", "rmerge@4"} <= set(first)
