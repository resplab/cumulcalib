#' Cumulative calibration assessment
#'
#' This is the core function for performing cumulative calibration assessment
#'
#' @return an objective of class cumulcalib that can be printed or plotted
#  @seealso [stringi::stri_length()] which this function wraps.
#' @param y vector of binary responses
#' @param p vector of predicted probabilities.
#' @param method string with either BB (Brownian bridge test, default method), BM (Brownian motion test), BM2p (two-part BM test - experimental), BB1p (one-part BB test wit only the 'bridge' component). Multiple methods can be specified. The first one will be the 'main' method (e.g., when submitting the resulting object to plot()). Default is c("BB","BM")
#' @param ordered if TRUE, y and p are already ordered based on ascending values of p. This is to speed up simulations.
#' @param n_sim if >0, indicates a simulation-based test is requested for inference.
#' @examples
#' pi <- rbeta(1000,1,2)
#' Y <- rbinom(length(pi),1,pi)
#' res <- cumulcalib(Y, pi, method="BB")
#' summary(res)
#' plot(res)
#' @export
cumulcalib <- function(
  y,
  p,
  method = c("BB", "BM"),
  ordered = FALSE,
  n_sim = 0
) {
  if (!ordered) {
    #Order ascendingly based on p, if not already ordered
    o <- order(p)
    p <- p[o]
    y <- y[o]
  }

  n <- length(p)

  #The time component of the random walk
  T_ <- sum(p * (1 - p)) #Total 'time'
  if (T_ < 30) {
    warning(
      "Total obsered time (sum(pi*(1-pi))) is less than 30; the data might be too small for reliable inference."
    )
  }
  t <- cumsum(p * (1 - p)) / T_ #time values at each p

  #Scaled cumulative sum (divided by n)
  C <- cumsum(y - p) / n #Scaled partial sum of prediction errors
  C_n <- C[n] #Mean calibration error
  C_star <- max(abs(C)) #Macimum distance

  #This is the S process as described in the paper. The function returns all these metrics regardless of which method is used. But inference is only done per specified method(s)
  scale <- n / sqrt(T_)
  S <- scale * C

  out <- inference(t, S, method)
  out$T <- T_
  out$C_n <- C_n
  out$C_star <- C_star
  out$data <- cbind(t = t, S = S, X = p, C = C) #Returns the generated random-walk
  class(out) <- c("cumulcalib")
  return(out)
}


#' Cumulative calibration assessment for individual treatment effect (ITE) predictions
#'
#' This is the core function for performing cumulative calibration assessment for
#' models that predict individual treatment effects (treatment benefit).
#'
#' @return an objective of class cumulcalib that can be printed or plotted
#' @param y vector of binary responses
#' @param h vector of predicted treatment benefits (the predicted reduction in outcome risk due to treatment).
#' @param a treatment indicator (1 if treated, 0 if control).
#' @param p optional vector of predicted baseline risks (the risk without treatment). If omitted (NULL), the marginal test is performed using observed event rates in the treated and control groups; if supplied, the conditional test is performed. Default is NULL.
#' @param method string with either BB (Brownian bridge test, default method), BM (Brownian motion test), BM2p (two-part BM test - experimental), BB1p (one-part BB test wit only the 'bridge' component). Multiple methods can be specified. The first one will be the 'main' method (e.g., when submitting the resulting object to plot()). Default is c("BB","BM")
#' @param ordered if TRUE, the data are already ordered based on ascending values of h. This is to speed up simulations.
#' @param n_sim if >0, indicates a simulation-based test is requested for inference.
#' @param aux if TRUE, auxiliary quantities (used internally and for diagnostics) are returned in the result. Default is FALSE.
#' @examples
#' p <- rbeta(1000, 1, 2)
#' a <- rbinom(length(p), 1, 0.5)
#' h <- runif(length(p), 0, 0.05)
#' Y <- rbinom(length(p), 1, pmax(0, pmin(1, p - a * h)))
#' res <- cumulcalibITE(Y, h = h, a = a)
#' summary(res)
#' plot(res)
#' @export
cumulcalibITE <- function(
  y,
  h,
  a,
  p = NULL,
  method = c("BB", "BM"),
  ordered = FALSE,
  n_sim = 0,
  aux = FALSE
) {
  if (!ordered) {
    #Order ascendingly based on Y, if not already ordered
    o <- order(h)
    h <- h[o]
    y <- y[o]
    a <- a[o]
    if (!is.null(p)) p <- p[o]
  }

  n <- length(y)
  k <- 1:n
  k1 <- cumsum(a)
  k0 <- k - k1
  Y1 <- cumsum(y)
  Y01 <- cumsum((1 - a) * y)
  Y11 <- cumsum(a * y)
  B <- k * (ifelse(k0 != 0, Y01 / k0, 0) - ifelse(k1 != 0, Y11 / k1, 0))

  if (!is.null(p)) {
    X_mu <- k *
      ((1 - a) *
        ifelse(k0 != 0, (y - p) / k0, 0) -
        a * ifelse(k1 != 0, (y - p + h) / k1, 0))
    C <- cumsum(X_mu) / n
    sigma2 <- k^2 *
      ((1 - a) *
        ifelse(k0 != 0, p * (1 - p) / k0^2, 0) +
        a * ifelse(k1 != 0, (p - h) * (1 - p + h) / k1^2, 0))
    s2 <- cumsum(sigma2)
    if (aux) {
      mu <- ifelse(
        a,
        c(0, Y01[-n]) /
          c(0, k0[-n]) -
          k * (c(0, Y11[-n]) + (p - h)) / k1 +
          (k - 1) * c(0, Y11[-n]) / (c(0, k1[-n])),
        (1 - a) *
          (k *
            (c(0, Y01[-n]) + p) /
            k0 -
            (k - 1) * c(0, Y01[-n]) / c(0, k0[-n]) -
            c(0, Y11[-n]) / c(0, k1[-n]))
      )
      mu[which(is.nan(mu))] <- 0
    }
  } else {
    C <- (B - cumsum(h)) / n
    s2 <- k^2 *
      (ifelse(k0 != 0, Y01 / k0 * (1 - Y01 / k0) / k0, 0) +
        ifelse(k1 != 0, Y11 / k1 * (1 - Y11 / k1) / k1, 0))
    if (aux) mu <- h
  }

  #The time component of the random walk
  T_ <- s2[n] #Total 'time'
  if (T_ < 30) {
    warning(
      "Total obsered time (sum(pi*(1-pi))) is less than 30; the data might be too small for reliable inference."
    )
  }
  t <- s2 / T_
  S <- C * n / sqrt(T_)

  #Inference part. We loop over the different methods requested by the user
  out <- inference(t, S, method)
  out$T <- T_
  out$C_n <- C[n] #Mean calibration error
  out$C_star <- max(abs(C)) #Maximum distance
  out$data <- cbind(t = t, S = S, B = B, X = h) #Returns the generated random-walk

  if (aux) {
    out$aux$mu <- mu
    out$aux$s2 <- s2
  }

  #Record which approach was used: marginal (no baseline risks) or conditional (predicted baseline risks supplied)
  out$approach <- if (is.null(p)) "marginal" else "conditional"

  class(out) <- c("cumulcalib", "cumulcalibITE")
  return(out)
}


inference <- function(t, S, method) {
  out <- list()
  methods <- list()
  n <- length(S)

  S_star <- max(abs(S))
  B_star <- max(abs(S - S[n] * t))
  S_n <- S[n]

  for (i in 1:length(method)) {
    mt <- method[i]

    if (mt %in% c('BM2p', 'BB')) {
      #Two-part BM and BB both generate component-specific p-values
      stat1 <- S[n]
      pval1 <- 2 * stats::pnorm(-abs(S[n]), 0, 1) #Two-sided z test for mean calibration

      if (mt == 'BB') {
        loc <- which.max(abs(S - S[n] * t)) #The bridge component of the BB test
        stat2 <- B_star
        pval2 <- 1 - pKolmogorov(stat2)
      } else {
        loc <- which.max(abs(S))
        stat2 <- max(abs(S))
        pval2 <- 1 - pMAD_BM_c(stat2, w1 = S[n])
      }

      fisher <- -2 * (log(pval1) + log(pval2)) #Fisher's method for combining p-values
      pval <- 1 - stats::pchisq(fisher, 4)

      methods[[mt]]$stat <- fisher
      methods[[mt]]$pval <- pval
      methods[[mt]]$stat_by_component <- c(mean = stat1, distance = stat2)
      methods[[mt]]$pval_by_component <- c(mean = pval1, distance = pval2)
      methods[[mt]]$loc <- loc
    }

    if (mt %in% c('BM', 'BB1p')) {
      #These are one-part methods
      if (mt == 'BB1p') {
        S2 <- S - S[n] * t
        loc <- which.max(abs(S2))
        stat <- max(abs(S2))
        pval <- 1 - pKolmogorov(stat)
        methods[[mt]]$stat <- stat
        methods[[mt]]$pval <- pval
        methods[[mt]]$loc <- loc
      } else {
        stat <- max(abs(S))
        loc <- which.max(abs(S))
        pval <- 1 - pMAD_BM(stat)
        methods[[mt]]$stat <- stat
        methods[[mt]]$pval <- pval
        methods[[mt]]$loc <- loc
      }
    }
  }

  out$method <- names(methods[1]) #The base method is the first requested one
  out$S_n <- S_n
  out$S_star <- S_star
  out$B_star <- B_star
  #Copy the first method results to root of the list
  for (nm in names(methods[[1]])) {
    out[nm] <- methods[[1]][nm]
  }
  out$by_method <- methods

  out
}


#' CDF of the distribution of the maximum absolute deviation of Brownian motion in \[0,1\] interval
#' @return a scalar value
#' @param q the quantity at which CDF will be evaluated. Currently accepts only a scalar
#' @param summands maximum number of terms to be evaluated in the infinite series (default=100)
#' @export
pMAD_BM <- function(q, summands = 100) {
  pi <- base::pi
  m <- 0:summands
  out <- sum((-1)^m / (2 * m + 1) * exp(-(2 * m + 1)^2 * pi^2 / 8 / q^2))

  return(4 / pi * out)
}


#' Quantile function of the distribution of the maximum absolute deviation of Brownian motion in \[0,1\] interval
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @export
qMAD_BM <- function(p) {
  x <- stats::uniroot(
    function(x) {
      pMAD_BM(x) - p
    },
    interval = c(0, 10)
  )
  unname(x$root)
}


#
# W2CDF <- function(q, n_terms=10)
# {
#   n <- 0:n_terms
#   A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
#   A2 <- erfc((4*n+1)/(2*sqrt(2*q)))
#
#   sqrt(2)*sum(A1*A2)
# }
#
#
# W2PDF <- function(x, n_terms=10)
# {
#   n <- 0:n_terms
#   A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
#   A2 <- (4*n+1)*exp(-((4*n+1)^2)/(8*x))
#
#   1/(2*sqrt(base:::pi*x^3))*sum(A1*A2)
# }

#' CDF of the Kolmogorov distribution
#' @return a scalar value
#' @param q the quantity at which CDF will be evaluated. Currently accepts only a scalar
#' @param summands maximum number of terms to be evaluated in the infinite series (default=ceiling(q*sqrt(72)+3/2))
#' @export
pKolmogorov <- function(q, summands = ceiling(q * sqrt(72) + 3 / 2)) {
  #if(!is.null(summands)) q <- q*(1+0.12/sqrt(summands)+0.11/summands)

  sqrt(2 * pi) *
    sapply(q, function(x) {
      if (x > 0) {
        sum(exp(-(2 * (1:summands) - 1)^2 * pi^2 / (8 * x^2))) / x
      } else {
        0
      }
    })
}


#' Quantile function of the Kolmogorov distribution
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @export
qKolmogorov <- function(p) {
  x <- stats::uniroot(
    function(x) {
      pKolmogorov(x) - p
    },
    interval = c(0, 10)
  )
  unname(x$root)
}


#' CDF of the distribution of the maximum absolute deviation of Brownian motion in \[0,1\] interval, conditional on its terminal value
#' @return a scalar value
#' @param q the quantity at which CDF will be evaluated. Currently accepts only a scalar
#' @param w1 the terminal value
#' @param method different infinite series to use (1,2,3)
#' @param exp_tolerance numerical tolerance as the stopping rule when evaluating the infinite sum (default -30 on the exponential scale)
#' @param summands number of terms to evaluate (default is ceiling(q * sqrt(72) + 3/2))
#' @export
pMAD_BM_c <- function(
  q,
  w1,
  method = 1,
  exp_tolerance = -30,
  summands = ceiling(q * sqrt(72) + 3 / 2)
) {
  if (q <= max(0, w1)) {
    return(0)
  }
  out <- c()
  if (1 %in% method) {
    A <- -2 * q * q
    B <- 2 * q * w1
    C <- -exp_tolerance
    D <- sqrt(B * B - 4 * A * C)

    n <- seq(round((-B + D) / A / 2, 0), round((-B - D) / A / 2, 0))
    un <- 2 * n * q

    terms <- un * w1 - un^2 / 2
    out <- c(out, sum((-1)^n * exp(terms)))
  }
  if (2 %in% method) {
    out <- c(
      out,
      sqrt(2 * pi) /
        q *
        sum(
          exp(w1^2 / 2 - (2 * (1:summands) - 1)^2 * pi^2 / (8 * q^2)) *
            (cos((2 * (1:summands) - 1) * pi * w1 / 2 / q))
        )
    )
  }
  if (3 %in% method) {
    #Wrong
    n <- -1000:1000
    x <- sum(
      stats::dnorm(w1 + 4 * n * q) - stats::dnorm(w1 + 4 * n * q + 2 * q)
    )
    out <- c(out, x)
  }

  out
}


#' Quantile function of the distribution of the maximum absolute deviation of Brownian motion in \[0,1\] interval, conditional on its terminal value
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @param w1 the terminal value
#' @export
qMAD_BM_c <- function(p, w1) {
  x <- stats::uniroot(
    function(x) {
      pMAD_BM_c(x, w1) - p
    },
    interval = c(0, 10)
  )
  unname(x$root)
}


# Human-readable label for a method code (internal helper)
.cumulcalib_method_label <- function(m) {
  switch(
    m,
    BM = "One-part Brownian motion (BM)",
    BM1p = "One-part Brownian motion (BM)",
    BM2p = "Two-part Brownian motion (BM2p)",
    BB1p = "One-part Brownian bridge (BB1p)",
    BB = "Two-part Brownian bridge (BB)",
    BB2p = "Two-part Brownian bridge (BB)",
    m
  )
}


# Direction and location of the maximum cumulative calibration error (C*).
# C* is the maximum absolute cumulative calibration error on the natural
# (risk/benefit) scale, and is the only reported quantity that admits a signed
# interpretation; the test statistics (S_n, S*, B*) are referred to null
# distributions of absolute deviations and are reported unsigned. The direction
# is derived on the fly from the returned process, using the same orientation
# (sign of the standardized cumulative sum at its peak) as plot.cumulcalib().
.cumulcalib_cstar_direction <- function(x) {
  S <- x$data[, 'S']
  loc <- which.max(abs(S))
  sgn <- sign(S[loc])
  is_ite <- inherits(x, "cumulcalibITE")
  quantity <- if (is_ite) "benefit" else "risk"
  phrase <- if (sgn > 0) {
    paste0("observed ", quantity, " > predicted")
  } else if (sgn < 0) {
    paste0("observed ", quantity, " < predicted")
  } else {
    "no net cumulative deviation"
  }
  list(
    loc = loc,
    time = unname(x$data[loc, 't']),
    pred = unname(x$data[loc, 'X']),
    sign = sgn,
    phrase = phrase,
    predictor_label = if (is_ite) "predicted ITE" else "predicted risk"
  )
}


# Shape of miscalibration around the maximum cumulative error (C*). At the C*
# location the cumulative process is at its extreme, so the local prediction
# error typically has opposite signs on the two sides. When the process returns
# meaningfully after the peak (an interior peak), this is a reversal/crossover;
# when the peak sits at an endpoint (monotone accumulation), the discrepancy is
# one-directional. `frac` is how large the post-peak return must be, relative to
# the peak, to count as a reversal; `min_frac` requires each side to hold a
# minimum share of the observations. Interpretation is only warranted when C* is
# itself notable; summary() gates the call on the standardized C* (S*).
.cumulcalib_crossover <- function(x, frac = 0.25, min_frac = 0.025) {
  S <- x$data[, 'S']
  n <- length(S)
  loc <- which.max(abs(S))
  left_net <- S[loc]
  right_net <- S[n] - S[loc]
  is_ite <- inherits(x, "cumulcalibITE")
  quantity <- if (is_ite) "benefit" else "risk"
  phr <- function(s) {
    if (s >= 0) paste0("observed ", quantity, " exceeds predicted")
    else paste0("observed ", quantity, " falls below predicted")
  }
  reverses <- (min(loc, n - loc) >= min_frac * n) &&
    (sign(right_net) != sign(left_net)) &&
    (abs(right_net) > frac * abs(left_net))
  list(
    reverses = reverses,
    pred     = unname(x$data[loc, 'X']),
    left     = phr(sign(left_net)),
    right    = phr(sign(right_net)),
    dominant = phr(sign(left_net)),
    predictor_label = if (is_ite) "predicted ITE" else "predicted risk"
  )
}


#' Concise printing of a cumulcalib object
#'
#' Prints a one-screen overview of a calibration assessment returned by
#' [cumulcalib()] or [cumulcalibITE()]. For a fuller breakdown use [summary()].
#'
#' @return The input object, invisibly.
#' @param x An object of class cumulcalib generated by cumulcalib() or cumulcalibITE()
#' @param ... Not used
#' @method print cumulcalib
#' @export
print.cumulcalib <- function(x, ...) {
  is_ite <- inherits(x, "cumulcalibITE")
  writeLines(if (is_ite) {
    "Cumulative calibration assessment (individualized treatment effects)"
  } else {
    "Cumulative calibration assessment (predicted risks)"
  })

  m <- x$method
  meta <- paste0(
    "Method: ",
    .cumulcalib_method_label(m),
    "   |   n = ",
    nrow(x$data)
  )
  if (is_ite && !is.null(x$approach)) {
    meta <- paste0("Approach: ", x$approach, "   |   ", meta)
  }
  writeLines(paste0("  ", meta))

  writeLines(sprintf("  Mean calibration error (C_n): %.4g", x$C_n))
  dir <- .cumulcalib_cstar_direction(x)
  writeLines(sprintf(
    "  Maximum cumulative calibration error (C*): %.4g  (%s, at %s = %.4g)",
    x$C_star,
    dir$phrase,
    dir$predictor_label,
    dir$pred
  ))

  mm <- x$by_method[[m]]
  if (m %in% c("BB", "BB2p", "BM2p")) {
    writeLines(sprintf(
      "  Mean component:   S_n = %.3g,  p = %.3g",
      x$S_n,
      mm$pval_by_component[["mean"]]
    ))
    if (m %in% c("BB", "BB2p")) {
      writeLines(sprintf(
        "  Bridge component: B* = %.3g,  p = %.3g",
        x$B_star,
        mm$pval_by_component[["distance"]]
      ))
    } else {
      writeLines(sprintf(
        "  Distance comp.:   S* = %.3g,  p = %.3g",
        x$S_star,
        mm$pval_by_component[["distance"]]
      ))
    }
    writeLines(sprintf(
      "  Unified p-value (moderate calibration): %.3g",
      mm$pval
    ))
  } else {
    writeLines(sprintf(
      "  Test statistic = %.3g,  p-value = %.3g",
      mm$stat,
      mm$pval
    ))
  }
  invisible(x)
}


#' Summarize a cumulcalib object
#'
#' Builds a structured summary of a calibration assessment for the selected
#' method. The returned object has class \code{summary.cumulcalib} and is
#' displayed by [print.summary.cumulcalib()].
#'
#' @return An object of class \code{summary.cumulcalib}: a list holding the key
#'   statistics for the selected method (mean calibration, distance/bridge
#'   statistic, p-values) together with the maximum cumulative calibration
#'   error C* and, since C* is the only quantity that admits a signed
#'   interpretation, its direction (\code{C_star_sign},
#'   \code{C_star_direction}) and the location of that maximum
#'   (\code{C_star_time}, \code{C_star_pred}). The test statistics are reported
#'   unsigned. When the standardized maximum deviation (S*) is at least
#'   \code{shape_threshold}, a \code{crossover} element describes the shape of
#'   the miscalibration around the C* location: whether the cumulative error
#'   reverses (an interior peak, with opposite directions on either side) or is
#'   one-directional (a monotone accumulation), together with the direction(s).
#'   Otherwise \code{crossover} is \code{NULL}.
#' @param object An object of class cumulcalib generated by cumulcalib() or cumulcalibITE()
#' @param method Which method to use. Options are BB (Brownian bridge test), BM (Brownian motion test), BB1p (1-part Brownian bridge test), and BM2p (2-part Brownian bridge test). If unspecified, uses the default method used in the cumulcalib() call
#' @param shape_threshold Minimum standardized maximum deviation (S*) for the shape of the miscalibration around C* to be described. Below this, the cumulative process is treated as unremarkable and no shape is reported (default 1.5).
#' @param ... Not used
#' @method summary cumulcalib
#' @export
summary.cumulcalib <- function(object, method = NULL, shape_threshold = 1.5, ...) {
  if (is.null(method)) {
    method <- object$method
  } else {
    if (length(method) > 1) {
      stop("Only one method can be requested.")
    }
    if (!(method %in% names(object$by_method))) {
      stop(paste(
        "The requested method has not been provided by the original cumulcalib() call. You asked for",
        method,
        "but the submitted object has method(s)",
        paste(names(object$by_method), collapse = ",")
      ))
    }
  }

  is_ite <- inherits(object, "cumulcalibITE")
  dir <- .cumulcalib_cstar_direction(object)
  # Describe the shape only when the deviation is notable (gate on S*).
  crossover <- if (object$S_star >= shape_threshold) {
    .cumulcalib_crossover(object)
  } else {
    NULL
  }

  structure(
    list(
      type = if (is_ite) "ITE" else "risk",
      approach = object$approach,
      method = method,
      n = nrow(object$data),
      C_n = object$C_n,
      C_star = object$C_star,
      C_star_sign = dir$sign,
      C_star_direction = dir$phrase,
      C_star_time = dir$time,
      C_star_pred = dir$pred,
      S_n = object$S_n,
      S_star = object$S_star,
      B_star = object$B_star,
      crossover = crossover,
      stats = object$by_method[[method]],
      predictor_label = dir$predictor_label
    ),
    class = "summary.cumulcalib"
  )
}


#' Print a cumulcalib summary
#'
#' @return The input object, invisibly.
#' @param x An object of class summary.cumulcalib generated by summary()
#' @param ... Not used
#' @method print summary.cumulcalib
#' @export
print.summary.cumulcalib <- function(x, ...) {
  writeLines(if (x$type == "ITE") {
    "Moderate calibration assessment of individualized treatment effects"
  } else {
    "Moderate calibration assessment of predicted risks"
  })
  if (x$type == "ITE" && !is.null(x$approach)) {
    writeLines(paste("Approach:", x$approach))
  }

  writeLines(paste("C_n (mean calibration error):", x$C_n))
  writeLines(paste0(
    "C* (maximum cumulative calibration error): ",
    x$C_star,
    "  (", x$C_star_direction, ")"
  ))
  writeLines(paste0(
    "  Location of maximum cumulative error: time = ",
    x$C_star_time,
    ", ",
    x$predictor_label,
    " = ",
    x$C_star_pred
  ))
  if (!is.null(x$crossover)) {
    cr <- x$crossover
    if (cr$reverses) {
      writeLines(paste0(
        "  Shape: the cumulative error reverses around ",
        x$predictor_label, " = ", signif(cr$pred, 3),
        " (", cr$left, " below this point; ", cr$right, " above it)"
      ))
    } else {
      writeLines(paste0(
        "  Shape: one-directional, no reversal (", cr$dominant,
        " across the range)"
      ))
    }
  }

  m <- x$method
  st <- x$stats
  writeLines(paste("Method:", .cumulcalib_method_label(m)))
  if (m %in% c("BM", "BM1p")) {
    writeLines(paste(
      "S* (test statistic for cumulative calibration error):",
      x$S_star
    ))
    writeLines(paste("p-value:", st$pval))
  } else if (m == "BM2p") {
    writeLines(paste("S_n (Z score for mean calibration error):", x$S_n))
    writeLines(paste(
      "S* (test statistic for cumulative calibration error):",
      x$S_star
    ))
    writeLines(paste0(
      "Component-wise p-values: mean calibration=",
      st$pval_by_component[["mean"]],
      " | Distance=",
      st$pval_by_component[["distance"]]
    ))
    writeLines(paste("Combined p-value (Fisher's method):", st$pval))
  } else if (m == "BB1p") {
    writeLines(paste(
      "B* (test statistic for maximum absolute bridged calibration error):",
      x$B_star
    ))
    writeLines(paste("Test statistic value:", st$stat))
    writeLines(paste("p-value:", st$pval))
  } else if (m %in% c("BB", "BB2p")) {
    writeLines(paste("S_n (Z score for mean calibration error):", x$S_n))
    writeLines(paste(
      "B* (test statistic for maximum absolute bridged calibration error):",
      x$B_star
    ))
    writeLines(paste0(
      "Component-wise p-values: mean calibration=",
      st$pval_by_component[["mean"]],
      " | Distance (bridged)=",
      st$pval_by_component[["distance"]]
    ))
    writeLines(paste("Combined p-value (Fisher's method):", st$pval))
  }
  invisible(x)
}


ate <- function(y, a) {
  A <- sum((1 - a) * y)
  B <- sum(1 - a)
  C <- sum(a * y)
  D <- sum(a)
  ate <- A / B - C / D
  v <- A / B * (1 - A / B) / B + C / D * (1 - C / D) / D

  c(ate = ate, var = v)
}
