"""Rank-one update formulae must match the R oracle AND dense recomputation."""

import numpy as np
from conftest import mat, vec

from catci import merging


def test_get_index_matches_oracle_masks(oracle):
    dx = oracle["rank_one_updates"]["dx"]
    dy = oracle["rank_one_updates"]["dy"]
    for cc in oracle["rank_one_updates"]["cases"]:
        i1 = merging.get_index(cc["dimension"], cc["j1"], dx, dy)
        i2 = merging.get_index(cc["dimension"], cc["j2"], dx, dy)
        np.testing.assert_array_equal(i1.astype(int), np.asarray(cc["index1"], dtype=int))
        np.testing.assert_array_equal(i2.astype(int), np.asarray(cc["index2"], dtype=int))


def test_updates_match_oracle(oracle, shared_TS):
    T, Sigma = shared_TS
    dx = oracle["rank_one_updates"]["dx"]
    dy = oracle["rank_one_updates"]["dy"]
    normsq0 = float(np.sum(T ** 2))
    tr0 = float(np.trace(Sigma))
    tr20 = float(np.sum(Sigma ** 2))

    for cc in oracle["rank_one_updates"]["cases"]:
        i1 = merging.get_index(cc["dimension"], cc["j1"], dx, dy)
        i2 = merging.get_index(cc["dimension"], cc["j2"], dx, dy)

        new_T = merging.update_T(T, i1, i2)
        new_Sigma = merging.update_Sigma(Sigma, i1, i2)
        np.testing.assert_allclose(new_T, vec(cc["new_T"]), atol=1e-9, rtol=1e-7)
        np.testing.assert_allclose(new_Sigma, mat(cc["new_Sigma"]), atol=1e-9, rtol=1e-7)

        normsq = merging.update_normsq(normsq0, T, i1, i2)
        tr = merging.update_tr(tr0, Sigma, i1, i2)
        tr2 = merging.update_tr2(tr20, Sigma, i1, i2)
        # matches R's fast update
        assert abs(normsq - cc["fast"]["normsq"]) < 1e-7
        assert abs(tr - cc["fast"]["tr"]) < 1e-7
        assert abs(tr2 - cc["fast"]["tr2"]) < 1e-7
        # and matches dense recomputation from the merged (T, Sigma)
        assert abs(normsq - float(np.sum(new_T ** 2))) < 1e-7
        assert abs(tr - float(np.trace(new_Sigma))) < 1e-7
        assert abs(tr2 - float(np.sum(new_Sigma ** 2))) < 1e-7
