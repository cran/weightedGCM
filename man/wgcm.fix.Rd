% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WGCM.R
\name{wgcm.fix}
\alias{wgcm.fix}
\title{Weighted Generalised Covariance Measure (WGCM) With Fixed Weight Functions
Conditional Independence Test}
\usage{
wgcm.fix(
  X,
  Y,
  Z,
  regr.meth,
  regr.pars = list(),
  weight.num,
  weight.meth = "sign",
  nsim = 499
)
}
\arguments{
\item{X}{A (n x d_X) numeric matrix with n observations of d_X variables.}

\item{Y}{A (n x d_Y) numeric matrix with n observations of d_Y variables.}

\item{Z}{A (n x d_Z) numeric matrix with n observations of d_Z variables.}

\item{regr.meth}{One of "gam" and "xgboost" indicating the regression method
used to estimate the conditional expectations E[X|Z] and E[Y|Z].}

\item{regr.pars}{Optional additional regression parameters according to
GeneralisedCovarianceMeasure::comp.resids().}

\item{weight.num}{Number k_0 of weight functions per dimension of Z to be
used additionally to the constant weight function w(z) = 1. The total number
of weight functions will be 1 + k_0 * d_Z. In case of max(d_X, d_Y) > 1, the
same 1 + k_0 * d_Z weight functions are used for every combination of the
components of X and Y.}

\item{weight.meth}{String indicating the method to choose the weight functions.
Currently, only "sign" is implemented.}

\item{nsim}{Number of samples used to calculate the p-value using simulation.}
}
\value{
A p-value for the null hypothesis of conditional independence of X
and Y given Z.
}
\description{
The Weighted Generalised Covariance Measure (WGCM) with Fixed Weight Functions
is a test for conditional independence. It is a generalisation of the
Generalised Covariance Measure implemented in the R package
GeneralisedCovarianceMeasure.
}
\examples{
set.seed(1)
n <- 200
Z <- rnorm(n)
X <- Z + 0.3*rnorm(n)
Y1 <- Z + 0.3*rnorm(n)
Y2 <- Z + 0.3*rnorm(n) + 0.3*X
Y3 <- Z + 0.3*rnorm(n) + 0.15*X^2
wgcm.fix(X, Y1, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
wgcm.fix(X, Y2, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
wgcm.fix(X, Y3, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")

}
\references{
Please cite the following papers.
Cyrill Scheidegger, Julia Hoerrmann, Peter Buehlmann:
"The Weighted Generalised Covariance Measure"
\url{https://arxiv.org/abs/2111.04361}

Rajen D. Shah, Jonas Peters:
"The Hardness of Conditional Independence Testing and the Generalised
Covariance Measure"
\url{https://arxiv.org/abs/1804.07203}
}
