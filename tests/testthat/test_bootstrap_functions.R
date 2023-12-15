test_that("matrix_sqrt gives a valid matrix square root for degenerate matrices", {
  m1 <- tcrossprod(c(1,-1,-1,1))
  expect_equal(m1, crossprod(matrix_sqrt(m1)))
  m2 <- crossprod(matrix((1:8)/100, nrow = 2, ncol = 4))
  expect_equal(m2, crossprod(matrix_sqrt(m2)))
})
