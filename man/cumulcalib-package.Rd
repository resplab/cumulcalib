% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumulcalib-package.R
\docType{package}
\name{cumulcalib-package}
\alias{cumulcalib-package}
\title{cumulcalib: Cumulative Calibration Assessment for Prediction Models}
\description{
Tools for visualization of, and inference on, the calibration of prediction models on the cumulative domain. This provides a method for evaluating calibration of risk prediction models without having to group the data or use tuning parameters (e.g., loess bandwidth). This package implements the methodology described in Sadatsafavi and Patkau (2024) \doi{10.1002/sim.10138}. The core of the package is cumulcalib(), which takes in vectors of binary responses and predicted risks. The plot() and summary() methods are implemented for the results returned by cumulcalib().
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/resplab/cumulcalib}
}

}
\author{
\strong{Maintainer}: Mohsen Sadatsafavi \email{mohsen.sadatsafavi@ubc.ca} (\href{https://orcid.org/0000-0002-0419-7862}{ORCID})

}
\keyword{internal}
