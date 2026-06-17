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


#' Summarizes a cumulcalib object
#' @return None
#' @param object An object of class cumulcalib generated by cumulcalib()
#' @param ... Other parameters passed to summary()
#' @export
summary.cumulcalib <- function(object, ...) {
  UseMethod("summary", object)
}


#' Summarizes a cumulcalib object
#' @return None
#' @param object An object of class cumulcalib generated by cumulcalib()
#' @param method Which method to use. Options are BB (Brownian bridge test), BM (Brownian motion test), BB1p (1-part Brownian bridge test), and BM2p (2-part Brownian bridge test). If unspecified, returns the default method used in the cumulcalib() call
#' @param ... Other parameters passed to summary()
#' @method summary cumulcalib
summary.cumulcalib <- function(object, method = NULL, ...) {
  if (is.null(method)) {
    method <- object$method
  } else {
    if (length(method) > 1) {
      stop("Only one method can be equested.")
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

  n <- dim(object$data)[1]
  writeLines(paste("C_n (mean calibration error):", object$C_n))
  writeLines(paste(
    "C* (maximum absolute cumulative calibration error):",
    object$C_star
  ))

  if (method == "BM" || method == "BM1p") {
    writeLines("Method: One-part Brownian Motion (BM)")
    writeLines(paste(
      "S* (test statistic for cumulative calibration error):",
      object$S_star
    ))
    writeLines(paste("p-value:", object$by_method$BM$pval))
    writeLines(paste(
      "Location of maximum drift:",
      object$by_method$BM$loc,
      " | time value:",
      object$data[object$by_method$BM$loc, 't'],
      " | predictor value:",
      object$data[object$by_method$BM$loc, 'X']
    ))
  }
  if (method == "BM2p") {
    writeLines("Method: Two-part Brownian Motion (BM2p)")
    writeLines(paste("S_n (Z score for mean calibration error)", object$S_n))
    writeLines(paste(
      "S* (test statistic for cumulative calibration error):",
      object$S_star
    ))
    writeLines(paste0(
      "Component-wise p-values: mean calibration=",
      object$by_method$BM2p$pval_by_component[1],
      " | Distance=",
      object$by_method$BM2p$pval_by_component[2]
    ))
    writeLines(object$by_method$BM2p$pval_by_component)
    writeLines(paste(
      "Combined p-value (Fisher's method):",
      object$by_method$BM2p$pval
    ))
    writeLines(paste(
      "Location of maximum drift:",
      object$by_method$BM2p$loc,
      " | time value:",
      object$data[object$by_method$BM2p$loc, 't'],
      " | predictor value:",
      object$data[object$by_method$BM2p$loc, 'X']
    ))
  }
  if (method == "BB1p") {
    writeLines("Method: One-part Brownian Bridge (BB1p)")
    writeLines(paste(
      "B* (test statistic for maximum absolute bridged calibration error):",
      object$B_star
    ))
    writeLines(paste("Test statistic value:", object$by_method$BB1p$stat))
    writeLines(paste("p-value:", object$by_method$BB1p$pval))
    writeLines(paste(
      "Location of maximum drift:",
      object$by_method$BB1p$loc,
      " | time value:",
      object$data[object$by_method$BB1p$loc, 't'],
      " | predictor value:",
      object$data[object$by_method$BB1p$loc, 'X']
    ))
  }
  if (method == "BB" || method == "BB2p") {
    writeLines("Method: Two-part Brownian Bridge (BB)")
    writeLines(paste("S_n (Z score for mean calibration error)", object$S_n))
    writeLines(paste(
      "B* (test statistic for maximum absolute bridged calibration error):",
      object$B_star
    ))
    writeLines(paste0(
      "Component-wise p-values: mean calibration=",
      object$by_method$BB$pval_by_component[1],
      " | Distance (bridged)=",
      object$by_method$BB$pval_by_component[2]
    ))
    writeLines(paste(
      "Combined p-value (Fisher's method):",
      object$by_method$BB$pval
    ))
    writeLines(paste(
      "Location of maximum drift:",
      object$by_method$BB$loc,
      " | time value:",
      object$data[object$by_method$BB$loc, 't'],
      " | predictor value:",
      object$data[object$by_method$BB$loc, 'X']
    ))
  }
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
