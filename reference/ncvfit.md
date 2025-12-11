# Direct interface for nonconvex penalized regression (non-pathwise)

This function is intended for users who know exactly what they're doing
and want complete control over the fitting process: no standardization
is applied, no intercept is included, no path is fit. All of these
things are best practices for data analysis, so if you are choosing not
to do them, you are on your own – there is no guarantee that your
results will be meaningful. Some things in particular that you should
pay attention to:

- If your model has an intercept, it is up to you to (un)penalize it
  properly, typically by settings its corresponding element of
  `penalty.factor` to zero.

- You should provide initial values for the coefficients; in nonconvex
  optimization, initial values are very important in determining which
  local solution an algorithm converges to.

## Usage

``` r
ncvfit(
  X,
  y,
  init = rep(0, ncol(X)),
  r,
  xtx,
  penalty = c("MCP", "SCAD", "lasso"),
  gamma = switch(penalty, SCAD = 3.7, 3),
  alpha = 1,
  lambda,
  eps = 1e-05,
  max.iter = 1000,
  penalty.factor = rep(1, ncol(X)),
  warn = TRUE
)
```

## Arguments

- X:

  Design matrix; no intercept will be added, no standardization will
  occur (n x p matrix)

- y:

  Response vector (length n vector)

- init:

  Initial values for beta. Default: zero (length p vector)

- r:

  Residuals corresponding to `init`; these will be calculated if not
  supplied, but if they have already been calculated elsewhere, it is
  more efficient to pass them as an argument. WARNING: If you supply an
  incorrect value of `r`, the solution will be incorrect. (length n
  vector)

- xtx:

  X scales: the jth element should equal `crossprod(X[,j])/n`. These
  will be calculated if not supplied, but if they have already been
  calculated elsewhere, it is more efficient to pass them as an
  argument. In particular, if X is standardized, one should pass
  `xtx = rep(1, p)`. WARNING: If you supply an incorrect value of `xtx`,
  the solution will be incorrect. (length p vector)

- penalty:

  Penalty function to be applied, either "MCP" (default), "SCAD", or
  "lasso")

- gamma:

  Tuning parameter of the MCP/SCAD penalty, as in
  [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md);
  default is 3 for MCP and 3.7 for SCAD.

- alpha:

  Tuning paramter controlling the ridge component of penalty, as in
  [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md);
  default is 1 (meaning no ridge penalty)

- lambda:

  Regularization parameter value at which to estimate beta; must be
  scalar – for pathwise optimization, see
  [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md)

- eps:

  Convergence threshhold. The algorithm iterates until the RMSD for the
  change in linear predictors for each coefficient is less than eps.
  Default is 1e-4.

- max.iter:

  Maximum number of allowed iterations; if this number is reached,
  algorithm will terminate prior to convergence. Default: 1000.

- penalty.factor:

  Multiplicative factor for the penalty applied to each coefficient, as
  in
  [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md).
  In particular, note that if you include an intercept, you probably
  want to set its entry to zero here.

- warn:

  Return warning messages for failures to converge and model saturation?
  Default is TRUE.

## Value

A list containing:

- `beta`: The estimated regression coefficients

- `iter`: The number of iterations required to solve for \`beta

- `loss`: The loss (residual sum of squares) at convergence

- `resid`: The residuals at convergence

- `lambda`: See above

- `penalty`: See above

- `gamma`: See above

- `alpha`: See above

- `penalty.factor`: See above

- `n`: Sample size

## Details

At the moment, this function only works for least-squares loss
functions. Additional functionality for other loss functions (logistic,
Cox) is in development.

## Examples

``` r
data(Prostate)
X <- cbind(1, Prostate$X)
y <- Prostate$y
fit <- ncvfit(X, y, lambda=0.1, penalty.factor=c(0, rep(1, ncol(X)-1)))
fit$beta
#>                    lcavol      lweight          age         lbph          svi 
#>  2.268444208  0.677388754  0.000000000 -0.013317940  0.143711214  0.000000000 
#>          lcp      gleason        pgg45 
#>  0.000000000  0.000000000  0.005398707 
# Compare with:
coef(ncvreg(X, y), 0.1)
#> (Intercept)                  lcavol     lweight         age        lbph 
#>  -0.6973059   0.0000000   0.5387509   0.6382717   0.0000000   0.0000000 
#>         svi         lcp     gleason       pgg45 
#>   0.6102800   0.0000000   0.0000000   0.0000000 
# The unstandardized version makes little sense here, as it fails to account
# for differences in the scales of the predictors.
```
