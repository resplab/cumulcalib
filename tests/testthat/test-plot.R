# Tests for the plot.cumulcalib method. Drawing is sent to a null PDF device so
# nothing is written to disk.

test_that("plot works for every method without error", {
  set.seed(21)
  p <- rbeta(1500, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p, method = c("BB", "BM", "BB1p", "BM2p"))

  pdf(NULL)
  on.exit(dev.off())

  for (m in c("BB", "BM", "BB1p", "BM2p")) {
    expect_no_error(plot(res, method = m))
  }
  expect_no_error(plot(res))                         # default method
  expect_no_error(plot(res, stats_config = NULL))    # no statistic lines
  expect_no_error(plot(res, draw_polygon = FALSE, x2axis = FALSE, y2axis = FALSE))
})

test_that("plot works for cumulcalibITE objects", {
  set.seed(22)
  n <- 1500
  p <- rbeta(n, 1, 2)
  a <- rbinom(n, 1, 0.5)
  h <- runif(n, 0, 0.05)
  y <- rbinom(n, 1, pmax(0, pmin(1, p - a * h)))
  res <- cumulcalibITE(y, h = h, a = a)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(res))
})

test_that("plot rejects invalid or multiple methods", {
  set.seed(23)
  p <- rbeta(500, 1, 2)
  y <- rbinom(length(p), 1, p)
  res <- cumulcalib(y, p, method = "BB")

  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot(res, method = "NOPE"))
  expect_error(plot(res, method = c("BB", "BM")))
})
