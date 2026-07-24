# Tests for the distribution functions (Kolmogorov and maximum-absolute-deviation
# of Brownian motion). These are deterministic, so we test mathematical
# properties: valid CDF range, monotonicity, limiting behaviour, and that each
# quantile function inverts its CDF.

test_that("pKolmogorov is a valid, monotonic CDF with correct limits", {
  expect_gte(pKolmogorov(0.5), 0)
  expect_lte(pKolmogorov(3), 1)

  q <- seq(0.3, 3, by = 0.1)
  vals <- vapply(q, pKolmogorov, numeric(1))
  expect_true(all(diff(vals) >= 0))      # non-decreasing
  expect_true(all(vals >= 0 & vals <= 1))

  expect_lt(pKolmogorov(0.3), 0.05)      # small q -> near 0
  expect_gt(pKolmogorov(3), 0.99)        # large q -> near 1
})

test_that("qKolmogorov inverts pKolmogorov", {
  for (p in c(0.5, 0.8, 0.9, 0.95)) {
    q <- qKolmogorov(p)
    expect_equal(pKolmogorov(q), p, tolerance = 1e-4)
  }
})

test_that("pMAD_BM is a valid, monotonic CDF and qMAD_BM inverts it", {
  expect_gte(pMAD_BM(1), 0)
  expect_lte(pMAD_BM(6), 1)

  q <- seq(0.5, 4, by = 0.1)
  vals <- vapply(q, pMAD_BM, numeric(1))
  expect_true(all(diff(vals) >= 0))
  expect_true(all(vals >= 0 & vals <= 1))

  for (p in c(0.5, 0.8, 0.95)) {
    expect_equal(pMAD_BM(qMAD_BM(p)), p, tolerance = 1e-4)
  }
})

test_that("pMAD_BM_c is 0 below the terminal value and is inverted by qMAD_BM_c", {
  # CDF is 0 when the quantile cannot exceed the (non-negative) terminal value
  expect_equal(as.numeric(pMAD_BM_c(0.3, w1 = 0.5)), 0)
  expect_equal(as.numeric(pMAD_BM_c(-1, w1 = 0)), 0)

  val <- as.numeric(pMAD_BM_c(2, w1 = 0.5))
  expect_true(is.finite(val) && val >= 0 && val <= 1)

  for (p in c(0.5, 0.8, 0.95)) {
    q <- qMAD_BM_c(p, w1 = 0.5)
    expect_equal(as.numeric(pMAD_BM_c(q, w1 = 0.5)), p, tolerance = 1e-3)
  }
})
