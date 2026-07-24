# Tests for the summary.cumulcalib method (returns a structured object that is
# displayed by print.summary.cumulcalib).

test_that("summary returns a structured summary object for each method", {
  set.seed(31)
  p <- rbeta(800, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p, method = c("BB", "BM", "BB1p", "BM2p"))

  for (m in c("BB", "BM", "BB1p", "BM2p")) {
    s <- summary(res, method = m)
    expect_s3_class(s, "summary.cumulcalib")
    expect_equal(s$method, m)
    expect_equal(s$type, "risk")
    expect_output(print(s))                 # detailed display works
  }
})

test_that("summary of an ITE object records type and approach", {
  set.seed(32)
  n <- 800
  p <- rbeta(n, 1, 2)
  a <- rbinom(n, 1, 0.5)
  h <- runif(n, 0, 0.05)
  y <- rbinom(n, 1, pmax(0, pmin(1, p - a * h)))

  s_marg <- summary(cumulcalibITE(y, h = h, a = a))
  s_cond <- summary(cumulcalibITE(y, h = h, a = a, p = p))
  expect_equal(s_marg$type, "ITE")
  expect_equal(s_marg$approach, "marginal")
  expect_equal(s_cond$approach, "conditional")
  expect_output(print(s_marg), "treatment effects")
})

test_that("summary rejects invalid or multiple methods", {
  set.seed(33)
  p <- rbeta(300, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p, method = "BB")

  expect_error(summary(res, method = "NOPE"))
  expect_error(summary(res, method = c("BB", "BM")))
})

test_that("summary exposes C* direction/location, and it is method-independent", {
  set.seed(34)
  n <- 3000
  p <- rbeta(n, 1, 2)
  a <- rbinom(n, 1, 0.5)
  h <- rep(0.10, n)               # over-predicts benefit
  y <- rbinom(n, 1, p)
  res <- cumulcalibITE(y, h = h, a = a, method = c("BB", "BM"))

  s_bb <- summary(res, method = "BB")
  s_bm <- summary(res, method = "BM")

  expect_true(all(c("C_star_sign", "C_star_direction", "C_star_time",
                    "C_star_pred") %in% names(s_bb)))
  # C* is a method-independent metric: direction must not depend on the method
  expect_equal(s_bb$C_star_sign, s_bm$C_star_sign)
  expect_equal(s_bb$C_star_direction, s_bm$C_star_direction)
  expect_equal(s_bb$C_star_sign, -1)
  expect_match(s_bb$C_star_direction, "observed benefit < predicted")
  # the old method-specific location field is gone
  expect_null(s_bb$loc_pred)
})

test_that("summary describes a reversing miscalibration and gates on S*", {
  set.seed(71)
  n <- 5000
  p <- runif(n)
  # under-predicts risk for low p, over-predicts for high p -> up-then-down process
  true <- pmin(pmax(p + 0.15 * (0.5 - p), 0.001), 0.999)
  y <- rbinom(n, 1, true)
  res <- cumulcalib(y, p, method = "BM")

  s <- summary(res)
  expect_false(is.null(s$crossover))
  expect_true(s$crossover$reverses)
  expect_match(s$crossover$left, "exceeds predicted")
  expect_match(s$crossover$right, "falls below predicted")
  # gating: an impossibly high threshold suppresses the shape description
  expect_null(summary(res, shape_threshold = Inf)$crossover)
})

test_that("summary reports one-directional shape for monotone miscalibration", {
  set.seed(72)
  n <- 5000
  p <- runif(n)
  true <- pmin(p + 0.10, 0.999)   # uniform under-prediction of risk
  y <- rbinom(n, 1, true)
  res <- cumulcalib(y, p, method = "BM")
  s <- summary(res, shape_threshold = 0)   # force a description
  expect_false(s$crossover$reverses)
})
