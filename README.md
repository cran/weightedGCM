
<!-- README.md is generated from README.Rmd. Please edit that file -->

# weightedGCM

<!-- badges: start -->
<!-- badges: end -->

The package weightedGCM contains two functions implementing the two
versions of the Weighted Generalised Covariance Measure (WGCM)
conditional independence test described in Scheidegger, Hoerrmann and
Buehlmann (2021) “The Weighted Generalised Covariance Measure”
\<arXiv:2111.04361>. It is a generalisation of the Generalised
Covariance Measure (GCM) implemented in the package
‘GeneralisedCovarianceMeasure’ by Jonas Peters and Rajen D. Shah based
on Shah and Peters (2020) “The Hardness of Conditional Independence
Testing and the Generalised Covariance Measure” \<arXiv:1804.07203>.

-   *wgcm.fix* calculates a p-value for the null hypothesis of
    conditional independence based on the WGCM using several fixed
    weight functions.

-   *wgcm.est* calculates a p-value for the null hypothesis of
    conditional indepencence based on the WGCM using a single estimated
    weight function.

## Installation

You can install the released version of weightedGCM from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("weightedGCM")
```

## Example

We generate some data to use the two conditional independence test.

``` r
library(weightedGCM)
## Generate data
set.seed(1)
n <- 200
Z <- rnorm(n)
X <- Z + 0.3*rnorm(n)

## Y1 _||_ X | Z
Y1 <- Z + 0.3*rnorm(n)
## Y2 not _||_ X | Z
Y2 <- Z + 0.3*rnorm(n) + 0.3*X
## Y3 not _||_ X | Z
Y3 <- Z + 0.3*rnorm(n) + 0.15*X^2

## Test for conditional independence using wgcm.fix()
wgcm.fix(X, Y1, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
#> [1] 0.868
wgcm.fix(X, Y2, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
#> [1] 0.004
wgcm.fix(X, Y3, Z, regr.meth = "gam", weight.num = 7, weight.meth = "sign")
#> [1] 0.004

## Test for conditional independence using wgcm.est()
wgcm.est(X, Y1, Z, beta = 0.3, regr.meth = "gam")
#> [1] 0.2389116
wgcm.est(X, Y2, Z, beta = 0.3, regr.meth = "gam")
#> [1] 9.974068e-05
wgcm.est(X, Y3, Z, beta = 0.3, regr.meth = "gam")
#> [1] 0.004474768
```
