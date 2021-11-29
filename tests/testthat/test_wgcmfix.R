test_that("weight_matrix() gives the right result", {
  Z <- seq(0, 1, length.out = 24)
  W1 <- cbind(rep(1,24), c(rep(-1,12), rep(1,12)))
  W2 <- cbind(rep(1,24), c(rep(-1,8), rep(1, 16)), c(rep(-1,16), rep(1,8)))
  W3 <- cbind(rep(1,24), c(rep(-1, 6), rep(1,18)),
              c(rep(-1,12), rep(1, 12)), c(rep(-1, 18), rep(1, 6)))
  W5 <- cbind(rep(1, 24), c(rep(-1, 4), rep(1, 20)), c(rep(-1, 8), rep(1, 16)),
              c(rep(-1, 12), rep(1, 12)), c(rep(-1, 16), rep(1, 8)), c(rep(-1, 20), rep(1, 4)))

  expect_equal(unname(weight_matrix(as.matrix(Z), weight.num = 1, weight.meth = "sign")), W1)
  expect_equal(unname(weight_matrix(as.matrix(Z), weight.num = 2, weight.meth = "sign")), W2)
  expect_equal(unname(weight_matrix(as.matrix(Z), weight.num = 3, weight.meth = "sign")), W3)
  expect_equal(unname(weight_matrix(as.matrix(Z), weight.num = 5, weight.meth = "sign")), W5)

  ZZ <- cbind(seq(0, 1, length.out = 6), - seq(0, 1, length.out = 6))
  W <- cbind(rep(1, 6), c(rep(-1, 2), rep(1, 4)), c(rep(-1, 4), rep(1, 2)), c(rep(1, 4), rep(-1, 2)), c(rep(1, 2), rep(-1,4)))
  expect_equal(unname(weight_matrix(as.matrix(ZZ), weight.num = 2, weight.meth= "sign")), W)
})


test_that("Z = NULL yields a warning", {
  #univariate
  set.seed(1)
  n <- 200
  X <- rnorm(n)
  Y <- X + 0.3*rnorm(n)
  expect_warning(wgcm.fix(X, Y, Z = NULL, regr.meth = "gam", weight.num = 7, weight.meth = "sign"),
                 "No Z specified. No weight functions can be calculated. Function simply tests for vanishing correlation between X and Y.")
  #multivariate
  set.seed(1)
  X <- cbind(rnorm(n), rnorm(n))
  Y <- cbind(X[,1] + 0.3*rnorm(n), rnorm(n))
  expect_warning(wgcm.fix(X, Y, Z = NULL, regr.meth = "gam", weight.num = 7, weight.meth = "sign"),
                 "No Z specified. No weight functions can be calculated. Function simply tests for vanishing correlation between components of X and Y.")

})

test_that("wgcm.fix works", {
  skip_on_cran()
  ##univariate X and Y
  set.seed(1)
  n <- 200
  Z <- rnorm(n)
  Zm <- cbind(Z, Z+rnorm(n), rnorm(n), rnorm(n), rnorm(n))
  X <- Z + 0.3*rnorm(n)
  Y <- Z + 0.3*rnorm(n)
  #weight.num = 7
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "xgboost", weight.num = 7, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "gam", weight.num = 7, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "xgboost", weight.num = 7, weight.meth = "sign"))

  #weight.num = 0
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "gam", weight.num = 0, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "xgboost", weight.num = 0, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "gam", weight.num = 0, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "xgboost", weight.num = 0, weight.meth = "sign"))

  ##multivariate X and Y
  X <- cbind(X, X+rnorm(n), rnorm(n))
  Y <- cbind(Y, rnorm(n))
  #weight.num = 7
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "xgboost", weight.num = 7, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "gam", weight.num = 7, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "xgboost", weight.num = 7, weight.meth = "sign"))
  #weight.num = 0
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "gam", weight.num = 0, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Z, regr.meth = "xgboost", weight.num = 0, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "gam", weight.num = 0, weight.meth = "sign"))
  expect_silent(wgcm.fix(X, Y, Zm, regr.meth = "xgboost", weight.num = 0, weight.meth = "sign"))
})



