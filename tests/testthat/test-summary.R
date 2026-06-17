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
