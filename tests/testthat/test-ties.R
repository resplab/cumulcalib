# Tests for the ties= argument of cumulcalib(), which handles observations
# with exactly tied predicted risk (p).

test_that("no ties: default behavior is unchanged and silent", {
  set.seed(101)
  p <- rbeta(500, 1, 2)              # continuous -> ties are (essentially) impossible
  y <- rbinom(length(p), 1, p)

  expect_no_message(res <- cumulcalib(y, p))
  expect_no_warning(cumulcalib(y, p))

  o <- order(p)
  expect_equal(res$C_n, mean((y - p)[o]), tolerance = 1e-10)
  expect_equal(res$C_star, max(abs(cumsum((y - p)[o]))) / length(p), tolerance = 1e-10)
})

test_that("ties = 'group' matches manual group-averaging of y", {
  p <- c(0.2, 0.5, 0.5, 0.5, 0.8, 0.8)
  y <- c(0,   1,   0,   0,   1,   0)
  n <- length(p)

  res <- suppressWarnings(suppressMessages(cumulcalib(y, p, ties = "group")))

  # Manual: replace y within each tied group of p by its group mean
  y_expected <- ave(y, p, FUN = mean)
  C_expected <- cumsum(y_expected - p) / n
  expect_equal(res$data[, "C"], C_expected, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(res$C_star, max(abs(C_expected)), tolerance = 1e-10)
})

test_that("ties = 'group' emits a message (not a warning) reporting the ties", {
  # Large enough that T_ = sum(p(1-p)) >= 30, so the unrelated small-sample
  # warning doesn't confound this check.
  set.seed(21)
  p <- c(rep(0.3, 60), rep(0.7, 60), runif(80, 0, 1))
  y <- rbinom(length(p), 1, p)

  expect_message(cumulcalib(y, p, ties = "group"), "2 groups of tied")
  expect_no_warning(suppressMessages(cumulcalib(y, p, ties = "group")))
})

test_that("ties = 'ignore' emits a warning and reproduces the legacy (no-averaging) computation", {
  p <- c(0.2, 0.5, 0.5, 0.5, 0.8, 0.8)
  y <- c(0,   1,   0,   0,   1,   0)
  n <- length(p)

  expect_warning(
    res <- cumulcalib(y, p, ties = "ignore"),
    "ties = \"ignore\""
  )
  # Legacy behavior: plain stable sort by p, no averaging, no reordering
  o <- order(p)
  C_expected <- cumsum((y - p)[o]) / n
  expect_equal(res$data[, "C"], C_expected, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("ties = 'random' emits a warning and is reproducible under a fixed seed", {
  p <- c(0.2, 0.5, 0.5, 0.5, 0.5, 0.8, 0.8)
  y <- c(0,   1,   0,   0,   1,   1,   0)

  expect_warning(
    { set.seed(42); res1 <- cumulcalib(y, p, ties = "random") },
    "randomly reordered"
  )
  set.seed(42)
  res2 <- suppressWarnings(cumulcalib(y, p, ties = "random"))

  expect_equal(res1$data, res2$data)
  expect_equal(res1$C_star, res2$C_star)
})

test_that("invalid ties argument is rejected", {
  p <- c(0.2, 0.5, 0.5, 0.8)
  y <- c(0, 1, 0, 1)
  expect_error(cumulcalib(y, p, ties = "nope"))
})

test_that("ties = 'group' is invariant to the input row order of tied observations, unlike 'ignore'", {
  # Construct one large tie block (p = 0.5) whose y values are deliberately
  # arranged in a way that is *not* exchangeable with respect to input row
  # order: all 1's, then all 0's. Under ties = "ignore" (stable sort keeps
  # ties in input order), this manufactures a spurious excursion inside the
  # tie block that has nothing to do with p -- exactly the failure mode
  # observed with real (large, tied) predicted-risk data.
  set.seed(11)
  n_tie <- 200
  y_tie_spiked <- c(rep(1, n_tie / 2), rep(0, n_tie / 2))

  p_before <- sort(runif(200, 0, 0.4))
  y_before <- rbinom(200, 1, p_before)
  p_after <- sort(runif(200, 0.6, 1))
  y_after <- rbinom(200, 1, p_after)

  build <- function(y_tie) {
    list(p = c(p_before, rep(0.5, n_tie), p_after), y = c(y_before, y_tie, y_after))
  }

  d_spiked <- build(y_tie_spiked)
  d_shuffled <- build(sample(y_tie_spiked)) # same multiset, different row order

  # "ignore": sensitive to the arbitrary within-tie input row order
  res_ig_1 <- suppressWarnings(cumulcalib(d_spiked$y, d_spiked$p, ties = "ignore"))
  res_ig_2 <- suppressWarnings(cumulcalib(d_shuffled$y, d_shuffled$p, ties = "ignore"))
  expect_false(isTRUE(all.equal(res_ig_1$C_star, res_ig_2$C_star)))

  # "group": identical regardless of within-tie input row order
  res_grp_1 <- suppressMessages(cumulcalib(d_spiked$y, d_spiked$p, ties = "group"))
  res_grp_2 <- suppressMessages(cumulcalib(d_shuffled$y, d_shuffled$p, ties = "group"))
  expect_equal(res_grp_1$data, res_grp_2$data)
  expect_equal(res_grp_1$C_star, res_grp_2$C_star)
  expect_equal(res_grp_1$by_method$BB$pval, res_grp_2$by_method$BB$pval)
})
