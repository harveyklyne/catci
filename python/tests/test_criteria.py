"""approx_chi CDF and the depth-0 comparators must match the R oracle."""

import numpy as np
from conftest import mat, vec

from catci import criteria


def test_approx_chi_matches_oracle(oracle):
    for cc in oracle["approx_chi"]["cases"]:
        got = criteria.approx_chi_metric(cc["normsq"], cc["tr"], cc["tr2"])
        assert abs(got - cc["value"]) < 1e-9


def test_scalar_methods_match_oracle(oracle, shared_TS):
    T, Sigma = shared_TS
    sm = oracle["scalar_methods"]
    assert abs(criteria.euclid(T, Sigma) - sm["euclid"]) < 1e-9
    assert abs(criteria.max_abs(T, Sigma) - sm["max"]) < 1e-9
    assert abs(criteria.mgcm(T, Sigma) - sm["mGCM"]) < 1e-9
