# Model predictions based on a fitted ncvreg object.

Similar to other predict methods, this function returns predictions from
a fitted `ncvreg` object.

## Usage

``` r
# S3 method for class 'cv.ncvreg'
predict(
  object,
  X,
  type = c("link", "response", "class", "coefficients", "vars", "nvars"),
  which = object$min,
  ...
)

# S3 method for class 'cv.ncvreg'
coef(object, which = object$min, ...)

# S3 method for class 'cv.ncvsurv'
predict(
  object,
  X,
  type = c("link", "response", "survival", "median", "coefficients", "vars", "nvars"),
  which = object$min,
  ...
)

# S3 method for class 'ncvreg'
predict(
  object,
  X,
  type = c("link", "response", "class", "coefficients", "vars", "nvars"),
  lambda,
  which = 1:length(object$lambda),
  ...
)

# S3 method for class 'ncvreg'
coef(object, lambda, which = 1:length(object$lambda), drop = TRUE, ...)
```

## Arguments

- object:

  Fitted `ncvreg` model object.

- X:

  Matrix of values at which predictions are to be made. Not used for
  `type="coefficients"` or for some of the `type` settings in `predict`.

- type:

  Type of prediction:

  - `link` returns the linear predictors

  - `response` gives the fitted values

  - `class` returns the binomial outcome with the highest probability

  - `coefficients` returns the coefficients

  - `vars` returns a list containing the indices and names of the
    nonzero variables at each value of `lambda`

  - `nvars` returns the number of nonzero coefficients at each value of
    `lambda`.

- which:

  Indices of the penalty parameter `lambda` at which predictions are
  required. By default, all indices are returned. If `lambda` is
  specified, this will override `which`.

- ...:

  Not used.

- lambda:

  Values of the regularization parameter `lambda` at which predictions
  are requested. For values of `lambda` not in the sequence of fitted
  models, linear interpolation is used.

- drop:

  If coefficients for a single value of `lambda` are to be returned,
  reduce dimensions to a vector? Setting `drop=FALSE` returns a 1-column
  matrix.

## Value

The object returned depends on type.

## References

Breheny P and Huang J. (2011) Coordinate descent algorithms for
nonconvex penalized regression, with applications to biological feature
selection. *Annals of Applied Statistics*, **5**: 232-253.
[doi:10.1214/10-AOAS388](https://doi.org/10.1214/10-AOAS388)

## See also

[`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md)

## Author

Patrick Breheny

## Examples

``` r
data(Heart)

fit <- ncvreg(Heart$X, Heart$y, family = "binomial")
coef(fit, lambda = 0.05)
#>  (Intercept)          sbp      tobacco          ldl    adiposity      famhist 
#> -4.079298688  0.000000000  0.037271649  0.075045442  0.000000000  0.611522063 
#>        typea      obesity      alcohol          age 
#>  0.009798506  0.000000000  0.000000000  0.047479496 
head(predict(fit, Heart$X, type = "link", lambda = 0.05))
#> [1]  0.358554123 -0.217849509 -0.510057641  0.546336589  0.216182502
#> [6] -0.007063715
head(predict(fit, Heart$X, type = "response", lambda = 0.05))
#> [1] 0.5886904 0.4457520 0.3751800 0.6332852 0.5538361 0.4982341
head(predict(fit, Heart$X, type = "class", lambda = 0.05))
#> [1] 1 0 0 1 1 0
predict(fit, type = "vars", lambda = c(0.05, 0.01))
#> $`0.0500`
#> tobacco     ldl famhist   typea     age 
#>       2       3       5       6       9 
#> 
#> $`0.0100`
#>     sbp tobacco     ldl famhist   typea obesity     age 
#>       1       2       3       5       6       7       9 
#> 
predict(fit, type = "nvars", lambda = c(0.05, 0.01))
#> 0.0500 0.0100 
#>      5      7 
```
