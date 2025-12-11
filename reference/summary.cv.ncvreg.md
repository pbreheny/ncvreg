# Summarizing cross-validation-based inference

Summary method for `cv.ncvreg` objects

## Usage

``` r
# S3 method for class 'cv.ncvreg'
summary(object, include_fit = FALSE, ...)

# S3 method for class 'summary.cv.ncvreg'
print(x, digits, include_fit = FALSE, ...)
```

## Arguments

- object:

  A `cv.ncvreg` or `cv.ncvsurv` object.

- include_fit:

  In addition to summarizing the cross-validation object, also summarize
  the model with the lowest CV error? (default: FALSE)

- ...:

  Further arguments passed to or from other methods.

- x:

  A `summary.cv.ncvreg` object.

- digits:

  Number of digits past the decimal point to print out. Can be a vector
  specifying different display digits for each of the five non-integer
  printed values.

## Value

An object with S3 class `summary.cv.ncvreg`. The class has its own print
method and contains the following list elements:

- penalty:

  The penalty used by `ncvreg`.

- model:

  Either `"linear"` or `"logistic"`, depending on the `family` option in
  `ncvreg`.

- n:

  Number of instances

- p:

  Number of regression coefficients (not including the intercept).

- min:

  The index of `lambda` with the smallest cross-validation error.

- lambda:

  The sequence of `lambda` values used by `cv.ncvreg`.

- cve:

  Cross-validation error (deviance).

- r.squared:

  Proportion of variance explained by the model, as estimated by
  cross-validation. For models outside of linear regression, the
  Cox-Snell approach to defining R-squared is used.

- snr:

  Signal to noise ratio, as estimated by cross-validation.

- sigma:

  For linear regression models, the scale parameter estimate.

- pe:

  For logistic regression models, the prediction error
  (misclassification error).

## References

Breheny P and Huang J. (2011) Coordinate descent algorithms for
nonconvex penalized regression, with applications to biological feature
selection. *Annals of Applied Statistics*, **5**: 232-253.
[doi:10.1214/10-AOAS388](https://doi.org/10.1214/10-AOAS388)

## See also

[`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md),
[`cv.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.md),
[`plot.cv.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/plot.cv.ncvreg.md),
[`summary.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/summary.ncvreg.md)

## Author

Patrick Breheny

## Examples

``` r
# Linear regression --------------------------------------------------
data(Prostate)
cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
summary(cvfit)
#> MCP-penalized linear regression with n=97, p=8
#> At minimum cross-validation error (lambda=0.0209):
#> -------------------------------------------------
#>   Nonzero coefficients: 7
#>   Cross-validation error (deviance): 0.57
#>   R-squared: 0.57
#>   Signal-to-noise ratio: 1.33
#>   Scale estimate (sigma): 0.752

# Logistic regression ------------------------------------------------
data(Heart)
cvfit <- cv.ncvreg(Heart$X, Heart$y, family="binomial")
summary(cvfit)
#> MCP-penalized logistic regression with n=462, p=9
#> At minimum cross-validation error (lambda=0.0117):
#> -------------------------------------------------
#>   Nonzero coefficients: 7
#>   Cross-validation error (deviance): 1.06
#>   R-squared: 0.20
#>   Signal-to-noise ratio: 0.26
#>   Prediction error: 0.277

# Cox regression -----------------------------------------------------
data(Lung)
cvfit <- cv.ncvsurv(Lung$X, Lung$y)
summary(cvfit)
#> MCP-penalized Cox regression with n=137, p=8
#> At minimum cross-validation error (lambda=0.2080):
#> -------------------------------------------------
#>   Nonzero coefficients: 2
#>   Cross-validation error (deviance): 7.59
#>   R-squared: 0.25
#>   Signal-to-noise ratio: 0.34
```
