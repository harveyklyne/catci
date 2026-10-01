"""The occupation taxonomy: a relabelling only, binary, contiguous, divisive-searchable."""

from __future__ import annotations

import numpy as np
import pytest

import adult
import adult_semisynth as ss
import taxonomies
from catci.search import MergeSearch, SplitSearch
from catci.structure import Ordinal


@pytest.fixture(scope="module", params=[1000, 2000, 30000])
def recoded(request):
    data = ss.pool_rare(adult.load(), "Occupation", 5, request.param)
    new, tree = taxonomies.apply_taxonomy(data, "Occupation")
    return data, new, tree


def _names(d):
    return d.codes["Occupation"].map(dict(enumerate(d.variables["Occupation"].levels, start=1)))


def test_recoding_is_a_relabelling(recoded):
    data, new, _ = recoded
    assert (_names(data) == _names(new)).all()
    others = [c for c in data.codes.columns if c != "Occupation"]
    assert new.codes[others].equals(data.codes[others])
    assert sorted(new.variables["Occupation"].levels) == sorted(data.variables["Occupation"].levels)


def test_tree_is_binary_with_contiguous_groups(recoded):
    _, new, tree = recoded
    d = len(new.variables["Occupation"].levels)
    stack = [tree.tree]
    while stack:
        node = stack.pop()
        labels = sorted(node.root)
        assert labels == list(range(labels[0], labels[-1] + 1))
        assert len(node.children) in (0, 2)
        stack.extend(node.children)
    assert sorted(tree.tree.root) == list(range(1, d + 1))


def test_first_split_is_white_collar(recoded):
    _, new, tree = recoded
    levels = new.variables["Occupation"].levels
    (first, _second), = tree.coarsest_partitions(len(levels))
    assert {levels[k - 1] for k in first} == {
        "Exec-managerial", "Prof-specialty", "Sales", "Tech-support", "Adm-clerical"}


def test_divisive_and_merge_search_run_on_the_tree(recoded):
    _, new, tree = recoded
    dx, dy = len(new.variables["Occupation"].levels), 2
    rng = np.random.default_rng(0)
    A = rng.normal(size=(dx * dy, dx * dy))
    Sigma = A @ A.T / (dx * dy)
    T = rng.normal(size=(dx * dy, 20))
    full = SplitSearch().paths(T, Sigma, dx, dy, tree, Ordinal())
    np.testing.assert_array_equal(SplitSearch(3).paths(T, Sigma, dx, dy, tree, Ordinal()), full[:4])
    assert MergeSearch().paths(T, Sigma, dx, dy, tree, Ordinal()).shape == full.shape
