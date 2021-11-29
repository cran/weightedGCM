## function for multivariate wgcm.fix
wgcm.fix.mult <- function(X, Y, Z, regr.meth, regr.pars, weight.num,
                          weight.meth, nsim){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- NROW(X)
  dx <- NCOL(X)
  dy <- NCOL(Y)
  if(is.null(Z)){
    warning("No Z specified. No weight functions can be calculated. Function simply tests for vanishing correlation between components of X and Y.")
    eps.mat <- t(t(X)-colMeans(X))
    xi.mat <- t(t(Y)-colMeans(Y))
    W <- rep(1,n)
  } else {
    Z <- as.matrix(Z)
    calc.res.Z <- function(V){
      return(as.numeric(GeneralisedCovarianceMeasure::comp.resids(V, Z,
                                                                  regr.pars = regr.pars,
                                                                  regr.method = regr.meth)))
    }
    eps.mat <- apply(X,2, calc.res.Z)
    xi.mat <- apply(Y,2,calc.res.Z)
    W <- weight_matrix(Z, weight.num, weight.meth)
  }
  R <- NULL
  for (j in 1:dx){
    for (l in 1:dy){
      R.jl <- eps.mat[,j]*xi.mat[,l]
      R <- cbind(R, R.jl*W)
    }
  }
  R <- t(R)
  R.norm <- R / sqrt(rowMeans(R^2) - rowMeans(R)^2)
  T.stat <- sqrt(n) * max(abs(rowMeans(R.norm)))
  T.stat.sim <- apply(abs(R.norm %*% matrix(rnorm(n * nsim), n, nsim)),
                      2, max) / sqrt(n)
  p.value <- (sum(T.stat.sim >= T.stat) + 1) / (nsim + 1)
  return(p.value)
}

## function for multivariate wgcm.est
wgcm.est.mult <- function(X, Y, Z, beta, regr.meth, regr.pars, nsim){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- NROW(X)
  dx <- NCOL(X)
  dy <- NCOL(Y)
  if(is.null(Z)){
    warning("No Z specified. No weight functions can be estimated. Function simply tests for vanishing correlation between components of X and Y.")
    eps.mat <- t(t(X)-colMeans(X))
    xi.mat <- t(t(Y)-colMeans(Y))
    R <- NULL
    for (j in 1:dx){
      for (l in 1:dy){
        R.jl <- eps.mat[,j]*xi.mat[,l]
        R <- cbind(R, R.jl)
      }
    }
    R <- t(R)
    R.norm <- R / sqrt(rowMeans(R^2) - rowMeans(R)^2)
    T.stat <- sqrt(n) * max(abs(rowMeans(R.norm)))
    T.stat.sim <- apply(abs(R.norm %*% matrix(rnorm(n * nsim), n, nsim)),
                        2, max) / sqrt(n)
    p.value <- (sum(T.stat.sim >= T.stat) + 1) / (nsim + 1)
  } else {
    Z <- as.matrix(Z)
    ind.train <- sample(1:n, ceiling(beta*n))
    Xtrain <- X[ind.train,]
    Xtest <- X[-ind.train,]
    Ytrain <- Y[ind.train,]
    Ytest <- Y[-ind.train,]
    Ztrain <- Z[ind.train,]
    Ztest <- Z[-ind.train,]
    calc.res.Z <- function(V){
      return(as.numeric(GeneralisedCovarianceMeasure::comp.resids(V, Ztest,
                                                                  regr.pars = regr.pars,
                                                                  regr.method = regr.meth)))
    }
    eps.mat <- apply(Xtest, 2, calc.res.Z)
    xi.mat <- apply(Ytest, 2, calc.res.Z)
    R <- NULL
    for (j in 1:dx){
      for (l in 1:dy){
        W <- predict_weight(Xtrain[,j], Ytrain[,l], Ztrain, Ztest, regr.meth, regr.pars)
        R.jl <- eps.mat[,j]*xi.mat[,l]
        R <- cbind(R, R.jl*W)
      }
    }
    R <- t(R)
    R.norm <- R/sqrt(rowMeans(R^2)-rowMeans(R)^2)
    ntest <- n-length(ind.train)
    #The estimated w-functions aim at making a positive test statistic. Hence
    #a one-sided test.
    T.stat <- sqrt(ntest) * max(rowMeans(R.norm))
    T.stat.sim <- apply(R.norm %*% matrix(rnorm(ntest * nsim), ntest, nsim),
                        2, max) / sqrt(ntest)
    p.value <- (sum(T.stat.sim >= T.stat) + 1)/(nsim + 1)
  }
  return(p.value)
}



## function to calculate weight matrix
weight_matrix <- function(Z, weight.num, weight.meth) {
  if (weight.meth == "sign") {
    n <- NROW(Z)
    dz <- NCOL(Z)
    W <- rep(1,n)
    if (weight.num >= 1) {
      d.probs <- (1:weight.num) / (weight.num + 1)
      for (i in 1:dz) {
        Zi <- Z[,i]
        a.vec <- quantile(Zi, d.probs, names=F)
        Wi <- outer(Zi, a.vec, signa)
        W <- cbind(W,Wi)
      }
    }
  } else {
    stop("Only method \"sign\" implemented yet to calculate weight function")
  }
  return(W)
}




## translated sign weight function
signa <- function(x,a){return(sign(x-a))}


## function to calculate a 1sided p-value for wgcm.fix, since we expect
## the test statistic to be positive under the alternative
wgcm.1d.1sided <- function(Xtest, Ytest, Ztest, W, regr.meth, regr.pars) {
  n <- NROW(Ztest)
  eps <- as.numeric(GeneralisedCovarianceMeasure::comp.resids(Xtest, Ztest,
                                                              regr.pars = regr.pars, regr.method = regr.meth))
  xi <- as.numeric(GeneralisedCovarianceMeasure::comp.resids(Ytest, Ztest,
                                                             regr.pars = regr.pars, regr.method = regr.meth))
  R <- eps*xi*W
  T.stat <- sqrt(n)*mean(R)/sqrt(mean(R^2)-mean(R)^2)
  p.value <- 1-pnorm(T.stat)
  return(p.value)
}

## function to estimate weight function for wgcm.est
predict_weight <- function(Xtrain, Ytrain, Ztrain, Ztest,
                           regr.meth, regr.pars) {
  eps <- as.numeric(GeneralisedCovarianceMeasure::comp.resids(Xtrain, Ztrain,
                                regr.pars = regr.pars, regr.method = regr.meth))
  xi <- as.numeric(GeneralisedCovarianceMeasure::comp.resids(Ytrain, Ztrain,
                                regr.pars = regr.pars, regr.method = regr.meth))
  eps.xi.train <- eps * xi
  switch(regr.meth, "gam"={
    W <- sign(predict_weight_gam(eps.xi.train, Ztrain, Ztest, regr.pars))
  }, "xgboost"={
    W <- sign(predict_weight_xgboost(eps.xi.train, Ztrain, Ztest, regr.pars))
  }
  )
  return(W)
}



## function to estimate weight function for wgcm.est using gam
predict_weight_gam <- function(V, Z, Ztest, pars = list()) {
  mod <- train.gam(Z, V, pars = list())
  dz <- NCOL(Z)
  dpred <- data.frame(Ztest)
  names(dpred) <- paste("var", 2:(dz+1), sep = "")
  Vpred <- predict(mod$model, dpred)
  return(Vpred)
}

## function to estimate weight function for wgcm.est using xgboost
predict_weight_xgboost <- function(V, Z, Ztest, pars = list()) {
  mod <- train.xgboost1(Z, V, pars = list())
  dz <- NCOL(Z)
  dpred <- xgb.DMatrix(data = data.matrix(Ztest))
  Vpred <- predict(mod$model, dpred)
  return(Vpred)
}


