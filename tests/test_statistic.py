"""approx_chi CDF and the depth-0 comparators must match the R oracle."""

import numpy as np
from conftest import mat, vec

from catci import statistic


def test_approx_chi_matches_oracle(oracle):
    for cc in oracle["approx_chi"]["cases"]:
        got = statistic.approx_chi_statistic(cc["normsq"], cc["tr"], cc["tr2"])
        assert abs(got - cc["value"]) < 1e-9


def test_scalar_methods_match_oracle(oracle, shared_TS):
    T, Sigma = shared_TS
    sm = oracle["scalar_methods"]
    assert abs(statistic.euclid(T, Sigma) - sm["euclid"]) < 1e-9
    assert abs(statistic.max_abs(T, Sigma) - sm["max"]) < 1e-9
    assert abs(statistic.mgcm(T, Sigma) - sm["mGCM"]) < 1e-9


def test_gammainc_is_bitwise_chi2_cdf(oracle):
    """approx_chi_statistic dropped chi2.cdf for speed; the values must not move at all."""
    from scipy.stats import chi2

    for cc in oracle["approx_chi"]["cases"]:
        g = cc["tr2"] / cc["tr"]
        h = cc["tr"] ** 2 / cc["tr2"]
        want = float(chi2.cdf(cc["normsq"] / g, df=h))
        assert statistic.approx_chi_statistic(cc["normsq"], cc["tr"], cc["tr2"]) == want
