"""greedy_search must reproduce the R greedy_query path (values + partitions).

This is the differential test that would have caught finding #1: the ordinal,
tree, and greedy searches must each match their own R path, and tree must differ
from ordinal.
"""

import numpy as np
from conftest import mat, vec

from catci.search import greedy_search
from catci.structure import Ordinal, Saturated, Tree


def _structures(name, dx, dy):
    if name == "ordinal":
        return Ordinal(), Ordinal()
    if name == "greedy":
        return Saturated(), Saturated()
    if name == "tree":
        return Tree.binary(dx), Tree.binary(dy)
    raise ValueError(name)


def _fixture_partitions(case):
    return [{"x": [list(g) for g in p["x"]], "y": [list(g) for g in p["y"]]}
            for p in case["partitions"]]


def test_search_paths_match_oracle(oracle, shared_TS):
    T, Sigma = shared_TS
    dx = oracle["search_paths"]["dx"]
    dy = oracle["search_paths"]["dy"]

    for name, case in oracle["search_paths"]["cases"].items():
        xs, ys = _structures(name, dx, dy)
        res = greedy_search(T, Sigma, dx, dy, xs, ys)

        np.testing.assert_allclose(
            np.asarray(res.values), vec(case["values"]), atol=1e-9, rtol=1e-7,
            err_msg=f"statistic values differ for method '{name}'",
        )
        assert _fixture_partitions(case) == res.partitions, f"partitions differ for method '{name}'"


def test_tree_differs_from_ordinal(oracle, shared_TS):
    # finding #1: genuine tree search must visit a different path than ordinal.
    T, Sigma = shared_TS
    dx = oracle["search_paths"]["dx"]
    dy = oracle["search_paths"]["dy"]
    ordv = greedy_search(T, Sigma, dx, dy, Ordinal(), Ordinal()).values
    treev = greedy_search(T, Sigma, dx, dy, Tree.binary(dx), Tree.binary(dy)).values
    assert len(ordv) == len(treev)
    assert not np.allclose(ordv, treev)
