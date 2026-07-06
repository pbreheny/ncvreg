# Extract Log-Likelihood

Extract the log-likelihood of an `ncvreg` or `ncvsurv` object.

## Usage

``` r
# S3 method for class 'ncvreg'
logLik(object, REML = FALSE, ...)

# S3 method for class 'ncvsurv'
logLik(object, ...)
```

## Arguments

- object:

  An `ncvreg` or `ncvsurv` object, as obtained from
  [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md) or
  [`ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/ncvsurv.md)

- REML:

  As in [`logLik.lm()`](https://rdrr.io/r/stats/logLik.html)

- ...:

  For S3 compatibility

## Value

An object with S3 class `logLik`. This is a vector of numbers (the
log-likelihood at each value of the regularization parameter) with the
following attributes:

- df:

  Degrees of freedom

- nobs:

  Number of observations

## See also

[`stats::logLik()`](https://rdrr.io/r/stats/logLik.html)
