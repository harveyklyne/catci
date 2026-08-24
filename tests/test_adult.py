"""Tests for the adult-income preprocessing.

The load-bearing test is :func:`test_fig1b_degrees_of_freedom_reproduce`. The
paper leaves the coding underspecified, but Fig. 1b prints six degrees of
freedom, and df is a deterministic function of the coding
(``df = (k-1)(r-1) * n_strata``). Reproducing all six pins the discretisation,
the level counts and the missing-value handling at once -- it is the only
checkable number the paper reports for this dataset.
"""

from __future__ import annotations

import numpy as np
import pytest

import adult


@pytest.fixture(scope="module")
def full():
    return adult.load()


def test_raw_file_is_the_uci_train_split():
    raw = adult.load_raw()
    assert len(raw) == 32561
    assert list(raw.columns) == adult.RAW_COLUMNS


def test_listwise_deletion_row_count(full):
    assert full.n == 30162
    assert adult.load(drop_missing=False).n == 32561


def test_eleven_analysis_variables(full):
    assert len(full.variables) == 11
    assert set(full.variables) == {
        "Age", "HoursPerWeek", "Education", "Income", "Workclass",
        "MaritalStatus", "Occupation", "Relationship", "Race", "Sex",
        "NativeCountry",
    }


def test_no_missing_category_survives(full):
    for name in ("Workclass", "Occupation", "NativeCountry"):
        assert "?" not in full.variables[name].levels


@pytest.mark.parametrize(
    "name,kind,n_levels",
    [
        ("Age", "ordinal", 7),           # paper: <21, 21-30, ..., 61-70, >70
        ("HoursPerWeek", "ordinal", 4),  # paper: <=20, 21-30, 30-40, >40
        ("Education", "ordinal", 16),
        ("Income", "ordinal", 2),        # binarised at $50K
        ("Relationship", "categorical", 6),
    ],
)
def test_declared_coding(full, name, kind, n_levels):
    assert full.kind(name) == kind
    assert full.n_levels(name) == n_levels


def test_codes_are_one_based_and_contiguous(full):
    for name in full.variables:
        col = full.codes[name].to_numpy()
        assert col.min() == 1
        assert col.max() == full.n_levels(name)
        assert set(np.unique(col)) == set(range(1, full.n_levels(name) + 1))


def test_fig1b_degrees_of_freedom_reproduce(full):
    """All six df of Fig. 1b, on a 1000-row subsample as the figure used.

    Workclass shows 6 levels and Occupation 13 rather than their full-data 7 and
    14: the rarest levels are absent from a sample this size, which is what makes
    the paper's numbers come out. Seed 6 is one such subsample.
    """
    check = adult.fig1b_df_check(full.subsample(1000, seed=6))
    assert check["match"].all(), check.to_string(index=False)
    assert (check["n_strata"] == 14).all()  # Age (7) x Sex (2)


def test_keeping_missing_as_a_category_cannot_reproduce_fig1b():
    """The evidence behind ``drop_missing=True`` being the default.

    With ``?`` retained, Workclass and Occupation carry an extra level, so no
    subsample can land on the 6 and 13 that Fig. 1b's df require.
    """
    kept = adult.load(drop_missing=False)
    for seed in range(10):
        check = adult.fig1b_df_check(kept.subsample(1000, seed=seed))
        assert not check["match"].all()


def test_subsample_recodes_to_contiguous_labels(full):
    sub = full.subsample(1000, seed=6)
    assert sub.n == 1000
    for name in sub.variables:
        col = sub.codes[name].to_numpy()
        assert set(np.unique(col)) == set(range(1, sub.n_levels(name) + 1))
        assert sub.n_levels(name) <= full.n_levels(name)


def test_subsample_is_deterministic(full):
    a = full.subsample(500, seed=3).codes
    b = full.subsample(500, seed=3).codes
    assert a.equals(b)
