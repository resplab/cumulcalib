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
cumulcalib <- function(y, p, method=c("BB","BM"), ordered=F, n_sim=0)
{
  out <- list()
  methods <- list()

  if(!ordered) #Order ascendingly based on p, if not already ordered
  {
    o <- order(p)
    p <- p[o]
    y <- y[o]
  }

  n <- length(p)

  #The time component of the random walk
  T_ <- sum(p*(1-p)) #Total 'time'
  if(T_<30) warning("Total obsered time (sum(pi*(1-pi))) is less than 30; the data might be too small for reliable inference.")
  t <- cumsum(p*(1-p))/T_ #time values at each p
  X <- p #Predicted values

  #Scaled cumulative sum (divided by n)
  C <- cumsum(y-p)/n  #Scaled partial sum of prediction errors
  C_n <- C[n] #Mean calibration error
  C_star <- max(abs(C))  #Macimum distance

  #This is the S process as described in the paper. The function returns all these metrics regardless of which method is used. But inference is only done per specified method(s)
  scale <- n/sqrt(T_)
  S <- scale*C
  S_star <- scale*C_star
  B_star <- max(abs(S-S[n]*t))
  S_n <- scale*C_n

  #Inference part. We loop over the different methods requested by the user
  for(i in 1:length(method))
  {
    mt <- method[i]

    if(mt %in% c('BM2p','BB')) #Two-part BM and BB both generate component-specific p-values
    {
      stat1 <- S[n]
      pval1 <- 2*stats::pnorm(-abs(S[n]),0,1) #Two-sided z test for mean calibration

      if(mt=='BB')
      {
        loc <- which.max(S-S[n]*t)  #The bridge component of the BB test
        stat2 <- B_star
        pval2 <- 1-pKolmogorov(stat2)
      }
      else
      {
        loc <- which.max(abs(S))
        stat2 <- max(abs(S))
        pval2 <- 1-pMAD_BM_c(stat2, w1=S[n])
      }

      fisher <- -2*(log(pval1)+log(pval2))  #Fisher's method for combining p-values
      pval <- 1-stats::pchisq(fisher,4)

      methods[[mt]]$stat <- fisher
      methods[[mt]]$pval <- pval
      methods[[mt]]$stat_by_component <- c(mean=stat1, distance=stat2)
      methods[[mt]]$pval_by_component <- c(mean=pval1, distance=pval2)
      methods[[mt]]$loc <- loc
    }

    if(mt %in% c('BM','BB1p'))  #These are one-part methods
    {
      if(mt=='BB1p')
      {
        S2 <- S-S[n]*t
        loc <- which.max(abs(S2))
        stat <- max(abs(S2))
        pval <- 1-pKolmogorov(stat)
        methods[[mt]]$stat <- stat
        methods[[mt]]$pval <- pval
        methods[[mt]]$loc <- loc
      }
      else
      {
        stat <- max(abs(S))
        loc <- which.max(abs(S))
        pval <- 1-pMAD_BM(stat)
        methods[[mt]]$stat <- stat
        methods[[mt]]$pval <- pval
        methods[[mt]]$loc <- loc
      }
    }
  }

  out$data <- cbind(X=X,t=t,C=C, S=S) #Returns the generated random-walk
  out$method <- names(methods[1]) #The base method is the first requested one

  out$scale <- scale
  out$C_n <- C_n
  out$S_n <- S_n
  out$C_star <- C_star
  out$S_star <- S_star
  out$B_star <- B_star

  #Copy the first method results to root of the list
  for(nm in names(methods[[1]]))
  {
    out[nm] <- methods[[1]][nm]
  }
  out$by_method <- methods

  class(out) <- "cumulcalib"
  return(out)
}





#' CDF of the distribution of the maximum absolute deviation of Brownian motion in \[0,1\] interval
#' @return a scalar value
#' @param q the quantity at which CDF will be evaluated. Currently accepts only a scalar
#' @param summands maximum number of terms to be evaluated in the infinite series (default=100)
#' @export
pMAD_BM <- function(q, summands=100)
{
  pi <- base::pi
  m<-0:summands
  out <- sum((-1)^m/(2*m+1)*exp(-(2*m+1)^2*pi^2/8/q^2))

  return(4/pi*out)
}


#' Quantile function of the distribution of the maximum absolute deviation of Brownian motion in \[0,1\] interval
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @export
qMAD_BM <- function(p)
{
  x <- stats::uniroot(function(x) {pMAD_BM(x)-p}, interval=c(0,10))
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
pKolmogorov <- function (q, summands=ceiling(q*sqrt(72)+3/2))
{
  #if(!is.null(summands)) q <- q*(1+0.12/sqrt(summands)+0.11/summands)

  sqrt(2 * pi) * sapply(q, function(x) {
    if (x > 0) {
      sum(exp(-(2 * (1:summands) - 1)^2 * pi^2/(8 * x^2)))/x
    }
    else {
      0
    }
  })
}


#' Quantile function of the Kolmogorov distribution
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @export
qKolmogorov <- function(p)
{
  x <- stats::uniroot(function(x) {pKolmogorov(x)-p}, interval=c(0,10))
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
pMAD_BM_c <- function(q, w1, method=1, exp_tolerance=-30, summands = ceiling(q * sqrt(72) + 3/2))
{
  if(q<=max(0,w1)) return(0)
  out <- c()
  if(1 %in% method)
  {
    A <- -2*q*q
    B <- 2*q*w1
    C <- -exp_tolerance
    D <- sqrt(B*B-4*A*C)

    n <- seq(round((-B+D)/A/2,0), round((-B-D)/A/2,0))
    un <- 2*n*q

    terms <- un*w1-un^2/2
    out <- c(out,sum((-1)^n*exp(terms)))
  }
  if(2 %in% method)
  {
    out <- c(out,sqrt(2*pi)/q*sum(exp(w1^2/2-(2*(1:summands)-1)^2*pi^2/(8*q^2))*(cos((2*(1:summands)-1)*pi*w1/2/q))))
  }
  if(3 %in% method)
  { #Wrong
    n <- -1000:1000
    x <- sum(stats::dnorm(w1+4*n*q)-stats::dnorm(w1+4*n*q+2*q))
    out <- c(out,x)
  }

  out
}




#' Quantile function of the distribution of the maximum absolute deviation of Brownian motion in \[0,1\] interval, conditional on its terminal value
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @param w1 the terminal value
#' @export
qMAD_BM_c <- function(p, w1)
{
  x <- stats::uniroot(function(x) {pMAD_BM_c(x,w1)-p}, interval=c(0,10))
  unname(x$root)
}






#' Generates cumulative calibration plots
#' @return None
#' @param x A cumulcalib object
#' @param ... Parameters to be passed to plot()
#' @export
plot.cumulcalib <- function(x,...)
{
  UseMethod("plot",x)
}

#' @return None
#' @param x An object of class cumulcalib generated by cumulcalib()
#' @param method Which method to use. Options are BB (Brownian bridge test), BM (Brownian motion test), BB1p (1-part Brownian bridge test), and BM2p (2-part Brownian bridge test). If unspecified, returns the default method used in the cumulcalib() call
#' @param draw_stat Should the statistic (terminal value an/or maximum drift, depending on method) be drawn? Default is TRUE
#' @param stat_col The color(s) to draw the stat. Default is c('blue','red'). For single-part tests (BM and BB1p) only the first element is used
#' @param draw_sig Whether significance lines should be drawn. Default is T. Colors will be the same as stat_col
#' @param sig_level If to draw significance lines, at what level? Default is c(0.95,0.95). For single-part tests (BM and BB1p) only the first element is used
#' @param x2axis If true, draws a second x-axis (on top) showing predicted risks
#' @param y2axis If true, draws a second y-axis (on right) showing scaled partial sums
#' @param ... Parameters to be passed to plot()
#' @method plot cumulcalib
plot.cumulcalib <- function(x, method=NULL, draw_stat=T, stat_col=c('blue','red'), draw_sig=T, sig_level=c(0.95,0.95), x2axis=T, y2axis=T, ...)
{
  args <- list(...)
  if(is.null(method))
  {
    method <- x$method
  }
  else
  {
    if(!("BM" %in% names(x$by_method)))
       stop("The requested method has not been requested in the original cumulcalib() call")
  }

  t_ <- x$data[,'t']
  X <- x$data[,'X']
  W <- x$data[,'S']

  n <- length(W)

  sign_p1 <- sign(W[n])
  sign_p2 <- sign(W[x$by_method[[method]]$loc])

  if(method=="BCI1p") #The only method that messes with S when drawing it.
  {
    W<-W-t_*W[n]
  }

  sig_p1 <- sig_p2 <- 0 #0 indicates do not draw signifcance lines
  if(draw_sig)
  {
    sig_p1 <- stats::qnorm(1-(1-sig_level[1])/2)

    if(method %in% c("BM"))
    {
      sig_p2 <- qMAD_BM(sig_level[2])
    }
    if(method %in% c("BM2p"))
    {
      sig_p2 <- qMAD_BM_c(sig_level[2],w1=unname(x$by_method$BM2p$stat_by_component['mean']))
    }
    if(method %in% c("BB1p","BB"))
    {
      sig_p2 <- qKolmogorov(sig_level[2])
    }
  }

  i <- match("xlim",names(args))
  if(is.na(i))
  {
    args$xlim<-c(-0.03,1.01)
  }
  i <- match("ylim",names(args))
  if(is.na(i))
  {
    args$ylim<- c(min(W,-sig_p1*(sign_p1==-1),-sig_p2*(sign_p2==-1),-1),max(W,sig_p1*(sign_p1==1),sig_p2*(sign_p2==1),1))
  }
  i <- match("type",names(args))
  if(is.na(i))
  {
    args$type<-'l'
  }
  i <- match("xaxt",names(args))
  if(is.na(i))
  {
    args$xaxt<-'n'
  }
  i <- match("yaxt",names(args))
  if(is.na(i))
  {
    args$yaxt<-'n'
  }
  i <- match("xlab",names(args))
  if(is.na(i))
  {
    args$xlab<-'Time'
  }
  i <- match("ylab",names(args))
  if(is.na(i))
  {
    args$ylab<-'Standardized cumulative sum'
  }

  args$x<-t_
  args$y<-W

  args$xaxs<-'i'

  par(mar=c(3,3,ifelse(x2axis,3,1),ifelse(y2axis,3,1)))

  do.call(plot, args)

  #Axes
  axis(1, line=0,  at=(0:10)/10, labels=(0:10)/10, padj=0)
  mtext(args$xlab, 1, line=2)

  axis(side=2, at=pretty(args$ylim), padj=0)
  mtext(args$ylab, 2, line=2)

  if(x2axis)
  {
    xs <- round(unlist(lapply((0:10)/10, function (z) {x$data[which.min((x$data[,'t']-z)^2),'X']})),3)
    axis(side=3, at=(0:10)/10, labels=xs, padj=0)
    mtext("Predicted values", 3, line=2) #expression(pi) for label will generate the symbol!
  }

  if(y2axis)
  {
    y2p <- pretty(args$ylim/x$scale)
    axis(side=4, at=y2p*x$scale, labels=y2p, padj=0)
    mtext('Scaled cumulative sum',4, line=2)
  }


  lines(c(0,1),c(0,0),col="grey")

  #Triangle
  {
    polygon(x=c(0,-0.03,-0.03), y=c(0,-1,1), col = 'black')
  }

  #P1 lines
  if(method %in% c("BM2p","BB"))
  {
    if(draw_stat) lines(c(1,1),c(0,W[n]), col=stat_col[1])
    if(draw_sig) lines(c(0,1),c(sign_p1*sig_p1,sign_p1*sig_p1),col=stat_col[1],lty=3)
  }

  #P2 lines
  loc <- x$by_method[[method]]$loc
  if(method %in% c('BB')) #If 2p bridge test then adjust the length of the red line and draw the bridge line
  {
    if(draw_stat) lines(c(0,1),c(0,W[n]),col="gray", lty=2)
    if(draw_stat) lines(c(t_[loc],t_[loc]),c(t_[loc]/t_[n]*W[n],W[loc]),col=stat_col[2])
    if(draw_sig) lines(c(0,1),c(sign_p2*sig_p2,sign_p2*sig_p2+W[n]),col=stat_col[2],lty=3)
  }
  else if(method %in% c('BCI1p'))
  {
    lines(c(0,1),c(0,W[n]),col="gray")
    if(draw_stat) lines(c(t_[loc],t_[loc]),c(0,W[loc]),col=stat_col[2])
    if(draw_sig) lines(c(0,1),c(sign_p2*sig_p2,sign_p2*sig_p2),col=stat_col[2],lty=3)
  }
  else #BM or BM2p
  {
    if(draw_stat) lines(c(t_[loc],t_[loc]),c(0,W[loc]),col=stat_col[2])
    if(draw_sig)  lines(c(0,1),c(sign_p2*sig_p2,sign_p2*sig_p2),col=stat_col[2],lty=3)
  }
}





#' Summarizes a cumulcalib object
#' @return None
#' @param object An object of class cumulcalib generated by cumulcalib()
#' @param ... Other parameters passed to summary()
#' @export
summary.cumulcalib <- function(object, ...)
{
  UseMethod("summary",object)
}



#' Summarizes a cumulcalib object
#' @return None
#' @param object An object of class cumulcalib generated by cumulcalib()
#' @param method Which method to use. Options are BB (Brownian bridge test), BM (Brownian motion test), BB1p (1-part Brownian bridge test), and BM2p (2-part Brownian bridge test). If unspecified, returns the default method used in the cumulcalib() call
#' @param ... Other parameters passed to summary()
#' @method summary cumulcalib
summary.cumulcalib <- function(object, method=NULL, ...)
{
  if(is.null(method))
  {
    method <- object$method
  }
  else
  {
    if(!("BM" %in% names(object$by_method)))
      stop("The requested method has not been requested in the original cumulcalib() call")
  }

  n <- dim(object$data)[1]
  writeLines(paste("C_n (mean calibration error):",object$C_n))
  writeLines(paste("C* (maximum absolute cumulative calibration error):",object$C_star))

  if(method=="BM" || method=="BM1p")
  {
    writeLines("Method: One-part Brownian Motion (BM)")
    writeLines(paste("S* (test statistic for cumulative calibration error):",object$S_star))
    writeLines(paste("p-value:",object$by_method$BM$pval))
    writeLines(paste("Location of maximum drift:",object$by_method$BM$loc,
                " | time value:", object$data[object$by_method$BM$loc,'t'],
                " | predictor value:", object$data[object$by_method$BM$loc,'X']))
  }
  if(method=="BM2p")
  {
    writeLines("Method: Two-part Brownian Motion (BM2p)")
    writeLines(paste("S_n (Z score for mean calibration error)",object$S_n))
    writeLines(paste("S* (test statistic for cumulative calibration error):",object$S_star))
    writeLines(paste0("Component-wise p-values: mean calibration=", object$by_method$BM2p$pval_by_component[1], " | Distance=", object$by_method$BM2p$pval_by_component[2]))
    writeLines(object$by_method$BM2p$pval_by_component)
    writeLines(paste("Combined p-value (Fisher's method):",object$by_method$BM2p$pval))
    writeLines(paste("Location of maximum drift:",object$by_method$BM2p$loc,
                " | time value:", object$data[object$by_method$BM2p$loc,'t'],
                " | predictor value:", object$data[object$by_method$BM2p$loc,'X']))
  }
  if(method=="BB1p")
  {
    writeLines("Method: One-part Brownian Bridge (BB1p)")
    writeLines(paste("B* (test statistic for maximum absolute bridged calibration error):",object$B_star))
    writeLines(paste("Test statistic value:",object$by_method$BB1p$stat))
    writeLines(paste("p-value:",object$by_method$BB1p$pval))
    writeLines(paste("Location of maximum drift:",object$by_method$BB1p$loc,
                " | time value:", object$data[object$by_method$BB1p$loc,'t'],
                " | predictor value:", object$data[object$by_method$BB1p$loc,'X']))
  }
  if(method=="BB" || method=="BB2p")
  {
    writeLines("Method: Two-part Brownian Bridge (BB)")
    writeLines(paste("S_n (Z score for mean calibration error)",object$S_n))
    writeLines(paste("B* (test statistic for maximum absolute bridged calibration error):",object$B_star))
    writeLines(paste0("Component-wise p-values: mean calibration=", object$by_method$BB$pval_by_component[1], " | Distance (bridged)=", object$by_method$BB$pval_by_component[2]))
    writeLines(paste("Combined p-value (Fisher's method):",object$by_method$BB$pval))
    writeLines(paste("Location of maximum drift:",object$by_method$BB$loc,
                " | time value:", object$data[object$by_method$BB$loc,'t'],
                " | predictor value:", object$data[object$by_method$BB$loc,'X']))
  }
}



