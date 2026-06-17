# Tests for cumulcalib() (calibration assessment of predicted risks).

test_that("cumulcalib returns a well-formed object", {
  set.seed(1)
  p <- rbeta(2000, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p)

  expect_s3_class(res, "cumulcalib")
  expect_true(all(
    c("C_n", "C_star", "S_n", "S_star", "B_star", "T", "data", "by_method", "method") %in%
      names(res)
  ))
  expect_equal(res$method, "BB")                       # default base method
  expect_setequal(names(res$by_method), c("BB", "BM"))
  expect_equal(ncol(res$data), 4L)
  expect_equal(nrow(res$data), length(p))
})

test_that("summary statistics match their mathematical definitions", {
  set.seed(2)
  p <- rbeta(500, 2, 3)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p)

  # C_n is the mean prediction error
  expect_equal(res$C_n, mean(y - p), tolerance = 1e-10)
  # S_n is the standardized terminal value
  expect_equal(res$S_n, sum(y - p) / sqrt(sum(p * (1 - p))), tolerance = 1e-8)
  # C_star is the maximum absolute cumulative error (after ordering by p), /n
  o <- order(p)
  expect_equal(res$C_star, max(abs(cumsum((y - p)[o]))) / length(p), tolerance = 1e-10)
})

test_that("ordered = TRUE matches ordered = FALSE on pre-sorted input", {
  set.seed(3)
  p <- sort(rbeta(300, 1, 1))
  y <- rbinom(length(p), 1, p)
  r1 <- cumulcalib(y, p, ordered = FALSE)
  r2 <- cumulcalib(y, p, ordered = TRUE)

  expect_equal(r1$S_n, r2$S_n)
  expect_equal(r1$by_method$BB$pval, r2$by_method$BB$pval)
  expect_equal(r1$by_method$BM$pval, r2$by_method$BM$pval)
})

test_that("each method returns a valid p-value", {
  set.seed(4)
  p <- rbeta(1000, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p, method = c("BB", "BM", "BB1p", "BM2p"))

  for (m in names(res$by_method)) {
    pv <- res$by_method[[m]]$pval
    expect_true(pv >= 0 && pv <= 1, info = paste("method", m))
  }
})

test_that("clear miscalibration yields a small p-value", {
  set.seed(6)
  n <- 3000
  truep <- rbeta(n, 1, 3)
  y <- rbinom(n, 1, truep)
  pbad <- pmin(0.99, truep * 3)             # systematically over-predict risk
  res <- cumulcalib(y, pbad)
  expect_lt(res$pval, 0.01)
})

test_that("small samples trigger the low-information warning", {
  set.seed(5)
  p <- rbeta(20, 1, 2)                       # sum(p(1-p)) is well under 30
  y <- rbinom(length(p), 1, p)
  expect_warning(cumulcalib(y, p), "less than 30")
})
