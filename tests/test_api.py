"""End-to-end smoke test of the public entry point (oracle-propensity path)."""

import numpy as np

from catci import catci_test
from catci.structure import Ordinal, Tree


def test_catci_test_runs_end_to_end():
    rng = np.random.default_rng(0)
    n, dx, dy = 300, 4, 4
    f = rng.uniform(0.2, 1.0, size=(n, dx)); f /= f.sum(1, keepdims=True)
    g = rng.uniform(0.2, 1.0, size=(n, dy)); g /= g.sum(1, keepdims=True)
    x = np.array([rng.choice(dx, p=f[i]) + 1 for i in range(n)])
    y = np.array([rng.choice(dy, p=g[i]) + 1 for i in range(n)])

    res = catci_test(x, y, Ordinal(), Ordinal(), f=f, g=g, n_boot=50, rng=rng)
    assert 0.0 <= res.p_value <= 1.0
    assert res.statistics[0] == res.statistics[0]  # not NaN
    # first partition is the fully-split one, last has 2 groups per dimension
    assert res.partitions[0]["x"] == [[1], [2], [3], [4]]
    assert len(res.partitions[-1]["x"]) == 2 and len(res.partitions[-1]["y"]) == 2


def test_requires_propensities_or_learner():
    x = np.array([1, 2, 3, 1]); y = np.array([1, 2, 3, 1])
    try:
        catci_test(x, y, Ordinal(), Ordinal())
    except ValueError:
        return
    raise AssertionError("expected ValueError when neither (f,g) nor (learner,z) given")
