# Tests for the ties= argument of cumulcalibITE(), conditional approach
# (p supplied). See tests/testthat/test-ties.R for the (structurally simpler)
# cumulcalib() risk-model case.

test_that("no ties: default 'average' matches 'ignore' and is silent (conditional)", {
  set.seed(201)
  n <- 500
  p <- rbeta(n, 1, 2)
  a <- rbinom(n, 1, 0.5)
  h <- runif(n, 0, 0.05) # continuous -> ties are (essentially) impossible
  y <- rbinom(n, 1, pmax(0, pmin(1, p - a * h)))

  expect_no_message(res <- cumulcalibITE(y, h = h, a = a, p = p))
  expect_no_warning(cumulcalibITE(y, h = h, a = a, p = p))

  res_ignore <- suppressWarnings(cumulcalibITE(y, h = h, a = a, p = p, ties = "ignore"))
  expect_equal(res$data, res_ignore$data)
  expect_equal(res$C_star, res_ignore$C_star)
})

test_that("ties = 'average' (conditional) matches the hand-derived group formula", {
  # A single tie group (h* = 0.5) with heterogeneous p, both arms, preceded by
  # some untied context so K0_end/K1_end/K_end are non-trivial.
  h_pre <- c(0.1, 0.2, 0.3)
  a_pre <- c(0, 1, 0)
  y_pre <- c(0, 1, 1)
  p_pre <- c(0.2, 0.4, 0.3)

  h_star <- 0.5
  a_tie <- c(0, 0, 1, 1)
  y_tie <- c(1, 0, 1, 0)
  p_tie <- c(0.30, 0.50, 0.45, 0.55) # heterogeneous within the tie

  h_post <- c(0.7, 0.8)
  a_post <- c(1, 0)
  y_post <- c(0, 1)
  p_post <- c(0.6, 0.35)

  h <- c(h_pre, rep(h_star, 4), h_post)
  a <- c(a_pre, a_tie, a_post)
  y <- c(y_pre, y_tie, y_post)
  p <- c(p_pre, p_tie, p_post)

  res <- suppressWarnings(suppressMessages(cumulcalibITE(y, h = h, a = a, p = p, ties = "average")))

  # Hand-derived closed form for the tie group's aggregate contribution
  # (see the design discussion this implements: evaluate the group at its own
  # end-state K_end/K0_end/K1_end, using each member's own p, no pooling).
  k <- seq_along(h)
  k1_all <- cumsum(a)
  k0_all <- k - k1_all
  end_idx <- length(h_pre) + 4 # last row of the tie group
  start_idx <- length(h_pre) # last row *before* the tie group
  K_end <- k[end_idx]
  K0_end <- k0_all[end_idx]
  K1_end <- k1_all[end_idx]

  ctrl <- a_tie == 0
  trt <- a_tie == 1
  dC_group <- K_end * (sum(y_tie[ctrl] - p_tie[ctrl]) / K0_end -
    sum(y_tie[trt] - p_tie[trt] + h_star) / K1_end)
  ds2_group <- K_end^2 * (sum(p_tie[ctrl] * (1 - p_tie[ctrl])) / K0_end^2 +
    sum((p_tie[trt] - h_star) * (1 - p_tie[trt] + h_star)) / K1_end^2)

  # Convert to the S/t scale the object actually reports (S = C*n/sqrt(T_),
  # t = s2/T_), using T_ as returned by the object itself.
  T_ <- res$T
  expected_dS <- dC_group / sqrt(T_)
  expected_dt <- ds2_group / T_

  observed_dS <- unname(res$data[end_idx, "S"] - res$data[start_idx, "S"])
  observed_dt <- unname(res$data[end_idx, "t"] - res$data[start_idx, "t"])

  expect_equal(observed_dS, expected_dS, tolerance = 1e-8)
  expect_equal(observed_dt, expected_dt, tolerance = 1e-8)
})

test_that("ties = 'average' (conditional) produces a straight line (equal steps) across a tie group", {
  set.seed(203)
  h_pre <- sort(runif(50, 0, 0.4))
  a_pre <- rbinom(50, 1, 0.5)
  p_pre <- runif(50, 0.1, 0.9)
  y_pre <- rbinom(50, 1, pmax(0, pmin(1, p_pre - a_pre * h_pre)))

  h_star <- 0.5
  m <- 10
  a_tie <- rbinom(m, 1, 0.5)
  p_tie <- runif(m, 0.1, 0.9)
  y_tie <- rbinom(m, 1, pmax(0, pmin(1, p_tie - a_tie * h_star)))

  h_post <- sort(runif(50, 0.6, 1))
  a_post <- rbinom(50, 1, 0.5)
  p_post <- runif(50, 0.1, 0.9)
  y_post <- rbinom(50, 1, pmax(0, pmin(1, p_post - a_post * h_post)))

  h <- c(h_pre, rep(h_star, m), h_post)
  a <- c(a_pre, a_tie, a_post)
  y <- c(y_pre, y_tie, y_post)
  p <- c(p_pre, p_tie, p_post)

  res <- suppressWarnings(suppressMessages(cumulcalibITE(y, h = h, a = a, p = p, ties = "average")))
  tie_rows <- (length(h_pre) + 1):(length(h_pre) + m)

  # location (S) and time (t) must both advance in exactly equal steps across
  # the tied block -- i.e. a straight line, hence no possible interior extremum
  dS <- diff(res$data[tie_rows, "S"])
  dt <- diff(res$data[tie_rows, "t"])
  expect_equal(dS, rep(dS[1], length(dS)), tolerance = 1e-8)
  expect_equal(dt, rep(dt[1], length(dt)), tolerance = 1e-8)
})

test_that("ties = 'average' (conditional) is invariant to input row order, unlike 'ignore'", {
  set.seed(11)
  n_tie <- 200
  h_star <- 0.03
  a_tie <- c(rep(0, 100), rep(1, 100))
  p_tie <- runif(n_tie, 0.2, 0.6) # heterogeneous within the tie
  # deliberately spiked: y=1 for the first half of each arm, y=0 for the second half
  y_tie_spiked <- c(rep(c(1, 0), each = 50), rep(c(1, 0), each = 50))

  h_before <- sort(runif(300, -0.05, h_star - 0.001))
  a_before <- rbinom(300, 1, 0.5)
  p_before <- runif(300, 0.1, 0.9)
  y_before <- rbinom(300, 1, pmax(0, pmin(1, p_before - a_before * 0.01)))

  h_after <- sort(runif(300, h_star + 0.001, 0.1))
  a_after <- rbinom(300, 1, 0.5)
  p_after <- runif(300, 0.1, 0.9)
  y_after <- rbinom(300, 1, pmax(0, pmin(1, p_after - a_after * 0.01)))

  build <- function(y_tie, a_tie, p_tie) {
    list(
      y = c(y_before, y_tie, y_after),
      h = c(h_before, rep(h_star, n_tie), h_after),
      a = c(a_before, a_tie, a_after),
      p = c(p_before, p_tie, p_after)
    )
  }

  ord <- sample(n_tie) # permute the tie block's rows as intact (y,a,p) triples
  d_spiked <- build(y_tie_spiked, a_tie, p_tie)
  d_shuffled <- build(y_tie_spiked[ord], a_tie[ord], p_tie[ord])

  res_ig_1 <- suppressWarnings(cumulcalibITE(d_spiked$y, h = d_spiked$h, a = d_spiked$a, p = d_spiked$p, ties = "ignore"))
  res_ig_2 <- suppressWarnings(cumulcalibITE(d_shuffled$y, h = d_shuffled$h, a = d_shuffled$a, p = d_shuffled$p, ties = "ignore"))
  expect_false(isTRUE(all.equal(res_ig_1$C_star, res_ig_2$C_star)))

  res_av_1 <- suppressMessages(cumulcalibITE(d_spiked$y, h = d_spiked$h, a = d_spiked$a, p = d_spiked$p, ties = "average"))
  res_av_2 <- suppressMessages(cumulcalibITE(d_shuffled$y, h = d_shuffled$h, a = d_shuffled$a, p = d_shuffled$p, ties = "average"))
  expect_equal(res_av_1$data, res_av_2$data)
  expect_equal(res_av_1$C_star, res_av_2$C_star)
  expect_equal(res_av_1$by_method$BB$pval, res_av_2$by_method$BB$pval)
})

test_that("ties = 'average' emits a message; 'ignore'/'random' emit warnings (conditional)", {
  set.seed(22)
  n <- 200
  p <- c(rep(0.3, 60), rep(0.7, 60), runif(80, 0, 1))
  a <- rbinom(n, 1, 0.5)
  h <- c(rep(0.01, 60), rep(0.02, 60), runif(80, 0, 0.05))
  y <- rbinom(n, 1, pmax(0, pmin(1, p - a * h)))

  expect_message(cumulcalibITE(y, h = h, a = a, p = p, ties = "average"), "2 groups of tied")
  expect_no_warning(suppressMessages(cumulcalibITE(y, h = h, a = a, p = p, ties = "average")))

  expect_warning(cumulcalibITE(y, h = h, a = a, p = p, ties = "ignore"), "ties = \"ignore\"")
  expect_warning(cumulcalibITE(y, h = h, a = a, p = p, ties = "random"), "randomly reordered")
})

test_that("ties = 'average' (marginal) matches 'ignore' and is silent when there are no ties", {
  set.seed(24)
  n <- 500
  a <- rbinom(n, 1, 0.5)
  h <- runif(n, 0, 0.05) # continuous -> ties are (essentially) impossible
  y <- rbinom(n, 1, 0.3)

  expect_no_message(res <- cumulcalibITE(y, h = h, a = a))
  expect_no_warning(cumulcalibITE(y, h = h, a = a))

  res_ignore <- suppressWarnings(cumulcalibITE(y, h = h, a = a, ties = "ignore"))
  expect_equal(res$data, res_ignore$data)
  expect_equal(res$C_star, res_ignore$C_star)
})

test_that("ties = 'average' (marginal) emits a message", {
  set.seed(23)
  n <- 200
  a <- rbinom(n, 1, 0.5)
  h <- c(rep(0.01, 60), rep(0.02, 60), runif(80, 0, 0.05))
  y <- rbinom(n, 1, 0.3)

  expect_message(cumulcalibITE(y, h = h, a = a, ties = "average"), "2 groups of tied")
  expect_no_warning(suppressMessages(cumulcalibITE(y, h = h, a = a, ties = "average")))
})

test_that("ties = 'average' (marginal) is invariant to input row order, unlike 'ignore'", {
  set.seed(25)
  n_tie <- 200
  h_star <- 0.03
  a_tie <- c(rep(0, 100), rep(1, 100))
  # deliberately spiked: y=1 for the first half of each arm, y=0 for the second half
  y_tie_spiked <- c(rep(c(1, 0), each = 50), rep(c(1, 0), each = 50))

  h_before <- sort(runif(300, -0.05, h_star - 0.001))
  a_before <- rbinom(300, 1, 0.5)
  y_before <- rbinom(300, 1, 0.3)

  h_after <- sort(runif(300, h_star + 0.001, 0.1))
  a_after <- rbinom(300, 1, 0.5)
  y_after <- rbinom(300, 1, 0.3)

  build <- function(y_tie, a_tie) {
    list(
      y = c(y_before, y_tie, y_after),
      h = c(h_before, rep(h_star, n_tie), h_after),
      a = c(a_before, a_tie, a_after)
    )
  }

  ord <- sample(n_tie) # permute the tie block's rows as intact (y,a) pairs
  d_spiked <- build(y_tie_spiked, a_tie)
  d_shuffled <- build(y_tie_spiked[ord], a_tie[ord])

  res_ig_1 <- suppressWarnings(cumulcalibITE(d_spiked$y, h = d_spiked$h, a = d_spiked$a, ties = "ignore"))
  res_ig_2 <- suppressWarnings(cumulcalibITE(d_shuffled$y, h = d_shuffled$h, a = d_shuffled$a, ties = "ignore"))
  expect_false(isTRUE(all.equal(res_ig_1$C_star, res_ig_2$C_star)))

  res_av_1 <- suppressMessages(cumulcalibITE(d_spiked$y, h = d_spiked$h, a = d_spiked$a, ties = "average"))
  res_av_2 <- suppressMessages(cumulcalibITE(d_shuffled$y, h = d_shuffled$h, a = d_shuffled$a, ties = "average"))
  expect_equal(res_av_1$data, res_av_2$data)
  expect_equal(res_av_1$C_star, res_av_2$C_star)
  expect_equal(res_av_1$by_method$BB$pval, res_av_2$by_method$BB$pval)
})

test_that("ties = 'random' (conditional) is reproducible under a fixed seed", {
  h <- c(0.05, 0.2, 0.2, 0.2, 0.2, 0.05)
  a <- c(0, 0, 1, 0, 1, 1)
  y <- c(0, 1, 0, 0, 1, 1)
  p <- c(0.3, 0.4, 0.5, 0.3, 0.55, 0.6) # kept so that p - h stays a valid probability

  set.seed(9)
  res1 <- suppressWarnings(cumulcalibITE(y, h = h, a = a, p = p, ties = "random"))
  set.seed(9)
  res2 <- suppressWarnings(cumulcalibITE(y, h = h, a = a, p = p, ties = "random"))

  expect_equal(res1$data, res2$data)
  expect_equal(res1$C_star, res2$C_star)
})
