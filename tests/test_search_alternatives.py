"""The TODO item 7 searches: beam, random path, and replaying a fixed path.

Deterministic identities are pinned against the R fixture through
``greedy_search`` (itself pinned by ``test_search.py``): a width-1 beam *is*
greedy, and replaying greedy's merges reproduces its values. A wide beam at
``d = 4`` is exhaustive, so it must match a brute-force max at every depth.
The calibration checks mirror ``test_properties.py``.
"""

import itertools

import numpy as np
import pytest

from catci.calibrate import adaptive_pvalue
from catci.search import beam_search, evaluate_path, greedy_search, random_merges
from catci.statistic import ApproxChi
from catci.structure import Ordinal, Saturated, Tree


def _structures(name, dx, dy):
    return {
        "ordinal": (Ordinal(), Ordinal()),
        "greedy": (Saturated(), Saturated()),
        "tree": (Tree.binary(dx), Tree.binary(dy)),
    }[name]


@pytest.mark.parametrize("name", ["ordinal", "greedy", "tree"])
def test_width_one_beam_is_greedy(oracle, shared_TS, name):
    T, Sigma = shared_TS
    dx, dy = oracle["search_paths"]["dx"], oracle["search_paths"]["dy"]
    xs, ys = _structures(name, dx, dy)
    g = greedy_search(T, Sigma, dx, dy, xs, ys)
    b = beam_search(T, Sigma, dx, dy, xs, ys, width=1)
    assert b.values == g.values  # bit-exact, not allclose
    assert b.partitions == g.partitions
    assert b.merges == g.merges


@pytest.mark.parametrize("name", ["ordinal", "greedy", "tree"])
def test_evaluate_path_replays_greedy(oracle, shared_TS, name):
    T, Sigma = shared_TS
    dx, dy = oracle["search_paths"]["dx"], oracle["search_paths"]["dy"]
    xs, ys = _structures(name, dx, dy)
    g = greedy_search(T, Sigma, dx, dy, xs, ys)
    e = evaluate_path(T, Sigma, dx, dy, g.merges)
    assert e.values == g.values
    assert e.partitions == g.partitions


def _set_partitions(labels):
    """All set partitions of ``labels`` (lists of lists)."""
    if not labels:
        yield []
        return
    first, rest = labels[0], labels[1:]
    for part in _set_partitions(rest):
        yield [[first]] + part
        for k in range(len(part)):
            yield part[:k] + [[first] + part[k]] + part[k + 1:]


def _value_of(T, Sigma, dx, dy, px, py):
    """Dense statistic value of the partition pair (px, py) -- no update formulae."""
    Mx = np.zeros((len(px), dx))
    for a, g in enumerate(px):
        Mx[a, [j - 1 for j in g]] = 1.0
    My = np.zeros((len(py), dy))
    for a, g in enumerate(py):
        My[a, [k - 1 for k in g]] = 1.0
    M = np.kron(My, Mx)  # X fastest, matching the index convention
    Tm, Sm = M @ T, M @ Sigma @ M.T
    stat = ApproxChi()
    return stat.value(stat.init(Tm, Sm))


def test_wide_beam_is_exhaustive_at_d4():
    rng = np.random.default_rng(3)
    dx = dy = 4
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    parts = list(_set_partitions(list(range(1, 5))))
    for _ in range(5):
        T = np.linalg.cholesky(Sigma) @ rng.standard_normal(p) + rng.standard_normal(p)
        b = beam_search(T, Sigma, dx, dy, Saturated(), Saturated(), width=10_000)
        g = greedy_search(T, Sigma, dx, dy, Saturated(), Saturated())
        for depth, value in enumerate(b.values):
            groups = dx + dy - depth
            best = max(
                _value_of(T, Sigma, dx, dy, px, py)
                for px, py in itertools.product(parts, parts)
                if len(px) + len(py) == groups and len(px) >= 2 and len(py) >= 2
            )
            assert np.isclose(value, best, rtol=1e-10, atol=1e-12)
            assert value >= g.values[depth] - 1e-12


def test_beam_dedups_partitions():
    # Two merge orders reach the same pair; a width-w beam must hold w distinct pairs.
    rng = np.random.default_rng(4)
    dx = dy = 5
    p = dx * dy
    Sigma = np.eye(p)
    T = rng.standard_normal(p)
    b = beam_search(T, Sigma, dx, dy, Saturated(), Saturated(), width=50)
    assert len(b.values) == dx + dy - 3


@pytest.mark.parametrize("name", ["ordinal", "greedy", "tree"])
def test_random_merges_are_permitted(name):
    rng = np.random.default_rng(5)
    dx, dy = 8, 4
    xs, ys = _structures(name, dx, dy)
    for _ in range(20):
        merges = random_merges(dx, dy, xs, ys, rng)
        assert len(merges) == dx + dy - 4
        partition = {"x": [[j] for j in range(1, dx + 1)], "y": [[k] for k in range(1, dy + 1)]}
        for dimension, i, j in merges:
            key, st = ("x", xs) if dimension == 1 else ("y", ys)
            assert (i, j) in st.permitted_merges(partition[key])
            partition[key][i - 1] += partition[key][j - 1]
            del partition[key][j - 1]
        assert len(partition["x"]) == 2 and len(partition["y"]) == 2


def test_random_path_values_are_dense_values():
    rng = np.random.default_rng(6)
    dx, dy = 5, 4
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    T = rng.standard_normal(p)
    res = evaluate_path(T, Sigma, dx, dy, random_merges(dx, dy, Saturated(), Saturated(), rng))
    for value, part in zip(res.values, res.partitions):
        assert np.isclose(value, _value_of(T, Sigma, dx, dy, part["x"], part["y"]))


# --------------------------------------------------------------------------- #
# calibration under a known Gaussian (loose, rep-limited: see test_properties)
# --------------------------------------------------------------------------- #
def _gaussian_pvalues(make_search, seed, reps=400):
    rng = np.random.default_rng(seed)
    dx = dy = 4
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    sqrtS = np.linalg.cholesky(Sigma)
    out = np.empty(reps)
    for r in range(reps):
        search = make_search(Sigma, dx, dy, rng)
        out[r] = adaptive_pvalue(
            sqrtS @ rng.standard_normal(p), Sigma, dx, dy, None, None,
            n_boot=100, rng=rng, search=search,
        )
    return out


def _check_uniform(pvals):
    assert abs(pvals.mean() - 0.5) < 0.06
    for alpha in (0.05, 0.10, 0.20):
        rate = float(np.mean(pvals < alpha))
        assert abs(rate - alpha) < 0.03, f"rejection rate {rate:.3f} at alpha={alpha}"


def test_beam_calibrates():
    def make(Sigma, dx, dy, rng):
        return lambda T: beam_search(T, Sigma, dx, dy, Saturated(), Saturated(), width=3).values

    _check_uniform(_gaussian_pvalues(make, seed=10))


def test_random_fixed_path_calibrates():
    # one path per replicate, shared by the observed and every bootstrap draw
    def make(Sigma, dx, dy, rng):
        merges = random_merges(dx, dy, Saturated(), Saturated(), rng)
        return lambda T: evaluate_path(T, Sigma, dx, dy, merges).values

    _check_uniform(_gaussian_pvalues(make, seed=11))


def test_random_fresh_path_calibrates():
    # a new path for every draw: valid because the paths are i.i.d. and independent of T
    def make(Sigma, dx, dy, rng):
        return lambda T: evaluate_path(
            T, Sigma, dx, dy, random_merges(dx, dy, Saturated(), Saturated(), rng)
        ).values

    _check_uniform(_gaussian_pvalues(make, seed=12))
