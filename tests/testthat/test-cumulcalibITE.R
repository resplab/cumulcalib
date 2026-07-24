# Tests for cumulcalibITE() (calibration assessment of predicted treatment benefit).

# Helper: generate a synthetic ITE validation set where the outcome risk is
# reduced by the predicted benefit h among the treated (a == 1).
make_ite <- function(n, seed) {
  set.seed(seed)
  p <- rbeta(n, 1, 2)
  a <- rbinom(n, 1, 0.5)
  h <- runif(n, 0, 0.05)
  y <- rbinom(n, 1, pmax(0, pmin(1, p - a * h)))
  list(y = y, a = a, h = h, p = p)
}

test_that("cumulcalibITE returns a well-formed object with the ITE subclass", {
  d <- make_ite(2000, 11)
  res <- cumulcalibITE(d$y, h = d$h, a = d$a)

  expect_s3_class(res, "cumulcalib")
  expect_s3_class(res, "cumulcalibITE")
  expect_true(all(c("S_n", "C_n", "C_star", "by_method", "data", "method") %in% names(res)))
})

test_that("marginal (no p) and conditional (with p) both yield valid p-values", {
  d <- make_ite(2000, 12)
  res_marginal    <- cumulcalibITE(d$y, h = d$h, a = d$a)
  res_conditional <- cumulcalibITE(d$y, h = d$h, a = d$a, p = d$p)

  for (res in list(res_marginal, res_conditional)) {
    pv <- as.numeric(res$pval_by_component["mean"])
    expect_true(pv >= 0 && pv <= 1)
  }
})

test_that("ordered = TRUE matches ordered = FALSE for ITE on pre-sorted input", {
  d <- make_ite(500, 13)
  o <- order(d$h)
  r1 <- cumulcalibITE(d$y[o], h = d$h[o], a = d$a[o], ordered = TRUE)
  r2 <- cumulcalibITE(d$y[o], h = d$h[o], a = d$a[o], ordered = FALSE)
  expect_equal(r1$S_n, r2$S_n)
  expect_equal(r1$C_n, r2$C_n)
})

test_that("aux = TRUE returns auxiliary quantities", {
  d <- make_ite(500, 14)
  res <- cumulcalibITE(d$y, h = d$h, a = d$a, aux = TRUE)
  expect_false(is.null(res$aux))
  expect_true(all(c("mu", "s2") %in% names(res$aux)))
})
