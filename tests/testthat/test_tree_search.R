test_that("tree search selects a different partition than ordinal search", {
  # Pins CODE_REVIEW.md finding #1: the "tree" method must actually use tree
  # search, not silently fall back to ordinal. A calibration test would NOT catch
  # this (both searches are valid tests and both calibrate); the diagnostic is
  # that tree and ordinal visit *different* merge paths on the same (T, Sigma).
  #
  # colsample_bylevel = 1, so the returned criterion vector and the final
  # partition are deterministic functions of (T, Sigma); the seed only fixes the
  # fixed input we build.
  set.seed(1)
  d <- 8L
  p <- d * d
  A <- matrix(rnorm(p * p), p, p)
  Sigma <- crossprod(A) / p + diag(p) # symmetric PSD
  T_vector <- as.numeric(rnorm(p))

  ord_vals <- query_lookup("ordinal")(T_vector, Sigma, dx = d, dy = d)
  tree_vals <- query_lookup("tree")(T_vector, Sigma, dx = d, dy = d)

  # Same length (13 at d = 8), but the criterion values must differ.
  expect_equal(length(ord_vals), length(tree_vals))
  expect_false(isTRUE(all.equal(ord_vals, tree_vals)))

  # And the two searches must end on genuinely different partitions.
  ord_full <- greedy_query(T_vector, Sigma, dx = d, dy = d,
                           metric = "approx_chi",
                           xsearch = "ordinal", ysearch = "ordinal",
                           colsample_bylevel = 1,
                           trees = list(NULL, NULL))
  tree_full <- greedy_query(T_vector, Sigma, dx = d, dy = d,
                            metric = "approx_chi",
                            xsearch = "tree", ysearch = "tree",
                            colsample_bylevel = 1,
                            trees = list(make_binary_tree(d), make_binary_tree(d)))
  ord_final <- ord_full$categories[[length(ord_full$categories)]]
  tree_final <- tree_full$categories[[length(tree_full$categories)]]
  expect_false(identical(ord_final, tree_final))
})

test_that("get_num_levels tree branch is guarded at d = 2", {
  # Pins the second half of finding #1: the tree branch must carry the same
  # (d > 2) guard as greedy/ordinal, so num_rows does not over-count once a
  # dimension collapses to two still-sibling groups.
  tree2 <- make_binary_tree(2L)
  cat2 <- as.list(seq(2L))
  expect_equal(get_num_levels("tree", 2L, tree2, cat2), 0)
})
