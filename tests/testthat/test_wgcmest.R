test_that("Z = NULL yields a warning", {
  #univariate
  set.seed(1)
  n <- 200
  X <- rnorm(n)
  Y <- X + 0.3*rnorm(n)
  expect_warning(wgcm.est(X, Y, Z = NULL, beta = 0.3, regr.meth = "gam", nsim = 499),
                 "No Z specified. No weight function can be estimated. Function simply tests for vanishing correlation between X and Y.")
  #multivariate
  set.seed(1)
  X <- cbind(rnorm(n), rnorm(n))
  Y <- cbind(X[,1] + 0.3*rnorm(n), rnorm(n))
  expect_warning(wgcm.est(X, Y, Z = NULL, beta = 0.3, regr.meth = "gam", nsim = 499),
                 "No Z specified. No weight functions can be estimated. Function simply tests for vanishing correlation between components of X and Y.")

})

test_that("wgcm.est works", {
  skip_on_cran()
  ##univariate X,Y
  set.seed(1)
  n <- 200
  Z <- rnorm(n)
  Zm <- cbind(Z, Z+rnorm(n), rnorm(n), rnorm(n), rnorm(n))
  X <- Z + 0.3*rnorm(n)
  Y <- Z + 0.3*rnorm(n)
  #univariate Z
  expect_silent(wgcm.est(X, Y, Z, beta = 0.3, regr.meth = "gam", nsim = 499))
  expect_silent(wgcm.est(X, Y, Z, beta = 0.3, regr.meth = "xgboost", nsim = 499))
  #multivariate Z
  expect_silent(wgcm.est(X, Y, Zm, beta = 0.3, regr.meth = "gam", nsim = 499))
  expect_silent(wgcm.est(X, Y, Zm, beta = 0.3, regr.meth = "xgboost", nsim = 499))

  ##multivariate X, Y
  X <- cbind(X, X+rnorm(n), rnorm(n))
  Y <- cbind(Y, rnorm(n))
  #univariate Z
  expect_silent(wgcm.est(X, Y, Z, beta = 0.3, regr.meth = "gam", nsim = 499))
  expect_silent(wgcm.est(X, Y, Z, beta = 0.3, regr.meth = "xgboost", nsim = 499))
  #multivariate Z
  expect_silent(wgcm.est(X, Y, Zm, beta = 0.3, regr.meth = "gam", nsim = 499))
  expect_silent(wgcm.est(X, Y, Zm, beta = 0.3, regr.meth = "xgboost", nsim = 499))
})
