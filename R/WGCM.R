#' Weighted Generalised Covariance Measure (WGCM) With Fixed Weight Functions
#' Conditional Independence Test
#'
#' The Weighted Generalised Covariance Measure (WGCM) with Fixed Weight Functions
#' is a test for conditional independence. It is a generalisation of the
#' Generalised Covariance Measure implemented in the R package
#' GeneralisedCovarianceMeasure.
#'
#' @param X A (n x d_X) numeric matrix with n observations of d_X variables.
#' @param Y A (n x d_Y) numeric matrix with n observations of d_Y variables.
#' @param Z A (n x d_Z) numeric matrix with n observations of d_Z variables.
#' @param regr.meth One of "gam" and "xgboost" indicating the regression method
#' used to estimate the conditional expectations E[X|Z] and E[Y|Z].
#' @param regr.pars Optional additional regression parameters according to
#' GeneralisedCovarianceMeasure::comp.resids().
#' @param weight.num Number k_0 of weight functions per dimension of Z to be
#' used additionally to the constant weight function w(z) = 1. The total number
#' of weight functions will be 1 + k_0 * d_Z. In case of max(d_X, d_Y) > 1, the
#' same 1 + k_0 * d_Z weight functions are used for every combination of the
#' components of X and Y.
#' @param weight.meth String indicating the method to choose the weight functions.
#' Currently, only "sign" is implemented.
#' @param nsim Number of samples used to calculate the p-value using simulation.
#'
#' @return A p-value for the null hypothesis of conditional independence of X
#' and Y given Z.
#'
#' @references Please cite the following papers.
#' Cyrill Scheidegger, Julia Hoerrmann, Peter Buehlmann:
#' "The Weighted Generalised Covariance Measure"
#' \url{https://arxiv.org/abs/2111.04361}
#'
#' Rajen D. Shah, Jonas Peters:
#' "The Hardness of Conditional Independence Testing and the Generalised
#' Covariance Measure"
#' \url{https://arxiv.org/abs/1804.07203}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' Z <- rnorm(n)
#' X <- Z + 0.3*rnorm(n)
#' Y1 <- Z + 0.3*rnorm(n)
#' Y2 <- Z + 0.3*rnorm(n) + 0.3*X
#' Y3 <- Z + 0.3*rnorm(n) + 0.15*X^2
#' wgcm.fix(X, Y1, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
#' wgcm.fix(X, Y2, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
#' wgcm.fix(X, Y3, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
#'
#' @export
#'
#' @importFrom stats pnorm quantile rnorm
wgcm.fix <- function(X, Y, Z, regr.meth, regr.pars = list(),
                     weight.num, weight.meth = "sign", nsim = 499) {
  if ((NCOL(X) != 1) || (NCOL(Y) != 1)) {
    p.value <- wgcm.fix.mult(X, Y, Z, regr.meth, regr.pars,
                             weight.num, weight.meth, nsim)
  } else{
    n <- NROW(X)
    if (is.null(Z)) {
      warning("No Z specified. No weight functions can be calculated. Function simply tests for vanishing correlation between X and Y.")
      weight.num <- 0
      eps <- X - mean(X)
      xi <- Y - mean(Y)
    } else {
      Z <- as.matrix(Z)
      eps <- as.numeric(GeneralisedCovarianceMeasure::comp.resids(X, Z,
                                                                  regr.pars = regr.pars, regr.method = regr.meth))
      xi <- as.numeric(GeneralisedCovarianceMeasure::comp.resids(Y, Z,
                                                                 regr.pars = regr.pars, regr.method = regr.meth))
    }
    if (weight.num == 0) {
      R <- eps * xi
      T.stat <- sqrt(n) * mean(R) / sqrt(mean(R^2) - mean(R)^2)
      p.value <- 2 * pnorm(-abs(T.stat))
    } else {
      W <- weight_matrix(Z, weight.num, weight.meth)
      R <- eps * xi * W
      R <- t(R)
      R.norm <- R / sqrt(rowMeans(R^2) - rowMeans(R)^2)
      T.stat <- sqrt(n) * max(abs(rowMeans(R.norm)))
      T.stat.sim <- apply(abs(R.norm %*% matrix(rnorm(n * nsim), n, nsim)),
                          2, max) / sqrt(n)
      p.value <- (sum(T.stat.sim >= T.stat) + 1) / (nsim + 1)
    }
  }
  return(p.value)
}

#' Weighted Generalised Covariance Measure (WGCM) With Estimated Weight Function
#' Conditional Independence Test
#'
#' The Weighted Generalised Covariance Measure (WGCM) with Estimated Weight Function
#' is a test for conditional independence. It is a generalisation of the
#' Generalised Covariance Measure implemented in the R package
#' GeneralisedCovarianceMeasure.
#'
#' @param X A (n x d_X) numeric matrix with n observations of d_X variables.
#' @param Y A (n x d_Y) numeric matrix with n observations of d_Y variables.
#' @param Z A (n x d_Z) numeric matrix with n observations of d_Z variables.
#' @param beta A real number between 0 and 1 indicating the fraction of the sample
#' used to estimate the weight function.
#' @param regr.meth One of "gam" and "xgboost" indicating the regression method
#' used to estimate the conditional expectations E[X|Z] and E[Y|Z] and the
#' weight function sign(E[(X-E[X|Z])(Y-E[Y|Z])|Z]).
#' @param regr.pars Optional additional regression parameters according to
#' GeneralisedCovarianceMeasure::comp.resids()
#' @param nsim Number of samples used to calculate the p-value using simulation.
#' Only used if max(d_X, d_Y) > 1.
#'
#' @return A p-value for the null hypothesis of conditional independence of X
#' and Y given Z.
#'
#' @references Please cite the following papers.
#' Cyrill Scheidegger, Julia Hoerrmann, Peter Buehlmann:
#' "The Weighted Generalised Covariance Measure"
#' \url{https://arxiv.org/abs/2111.04361}
#'
#' Rajen D. Shah, Jonas Peters:
#' "The Hardness of Conditional Independence Testing and the Generalised
#' Covariance Measure"
#' \url{https://arxiv.org/abs/1804.07203}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' Z <- rnorm(n)
#' X <- Z + 0.3*rnorm(n)
#' Y1 <- Z + 0.3*rnorm(n)
#' Y2 <- Z + 0.3*rnorm(n) + 0.3*X
#' Y3 <- Z + 0.3*rnorm(n) + 0.15*X^2
#' wgcm.est(X, Y1, Z, beta = 0.3, regr.meth = "gam")
#' wgcm.est(X, Y2, Z, beta = 0.3, regr.meth = "gam")
#' wgcm.est(X, Y3, Z, beta = 0.3, regr.meth = "gam")
#'
#' @export
#'
#' @import xgboost
#' @import mgcv
#' @importFrom methods show
#' @importFrom stats formula
#' @importFrom stats predict
wgcm.est <- function(X, Y, Z, beta = 0.3, regr.meth, regr.pars = list(), nsim=499) {
  if ((NCOL(X) != 1) || (NCOL(Y) != 1)) {
    p.value <- wgcm.est.mult(X, Y, Z, beta, regr.meth, regr.pars, nsim)
  } else {
    n <- NROW(X)
    if (is.null(Z)){
      warning("No Z specified. No weight function can be estimated. Function simply tests for vanishing correlation between X and Y.")
      eps <- X - mean(X)
      xi <- Y - mean(Y)
      R <- eps * xi
      T.stat <- sqrt(n) * mean(R) / sqrt(mean(R^2) - mean(R)^2)
      p.value <- 2 * pnorm(-abs(T.stat))
    } else {
      Z <- as.matrix(Z)
      ind.train <- sample(1:n, ceiling(beta * n))
      Ztrain <- Z[ind.train, ]
      Ztest <- Z[-ind.train, ]
      W <-  predict_weight(X[ind.train], Y[ind.train], Ztrain, Ztest,
                           regr.meth, regr.pars)
      p.value <- wgcm.1d.1sided(X[-ind.train],Y[-ind.train], Ztest, W,
                                regr.meth, regr.pars)
    }
  }
  return(p.value)
}
