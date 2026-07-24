# Tests for the print.cumulcalib method.

test_that("print is concise and does not dump the data matrix", {
  set.seed(41)
  p <- rbeta(2000, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p)

  out <- capture.output(printed <- print(res))
  # concise: a handful of lines, far fewer than the 2000-row data matrix
  expect_lt(length(out), 12)
  expect_true(any(grepl("predicted risks", out)))
  expect_true(any(grepl("Unified p-value", out)))
  # print returns its input invisibly
  expect_identical(printed, res)
})

test_that("print labels ITE objects and shows the approach", {
  set.seed(42)
  n <- 1500
  p <- rbeta(n, 1, 2)
  a <- rbinom(n, 1, 0.5)
  h <- runif(n, 0, 0.05)
  y <- rbinom(n, 1, pmax(0, pmin(1, p - a * h)))

  out_marg <- capture.output(print(cumulcalibITE(y, h = h, a = a)))
  out_cond <- capture.output(print(cumulcalibITE(y, h = h, a = a, p = p)))
  expect_true(any(grepl("treatment effects", out_marg)))
  expect_true(any(grepl("Approach: marginal", out_marg)))
  expect_true(any(grepl("Approach: conditional", out_cond)))
})

test_that("print works for one-part methods", {
  set.seed(43)
  p <- rbeta(800, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p, method = "BM")
  out <- capture.output(print(res))
  expect_true(any(grepl("Test statistic", out)))
})

test_that("print reports the direction of C* (worded tag) for ITE", {
  set.seed(51)
  n <- 3000
  p <- rbeta(n, 1, 2)
  a <- rbinom(n, 1, 0.5)
  h <- rep(0.10, n)              # model claims benefit of 0.10 ...
  y <- rbinom(n, 1, p)           # ... but treatment truly has ~no effect
  out <- capture.output(print(cumulcalibITE(y, h = h, a = a)))
  expect_true(any(grepl("Maximum cumulative calibration error", out)))
  expect_true(any(grepl("observed benefit < predicted", out)))
})

test_that("print reports the direction of C* for risk models", {
  set.seed(52)
  n <- 3000
  p <- rbeta(n, 1, 2)
  y <- rbinom(n, 1, p)
  out <- capture.output(print(cumulcalib(y = y, p = pmin(1, p + 0.15), method = "BM")))
  expect_true(any(grepl("observed risk < predicted", out)))
})
