# Extract residuals from a ncvreg or ncvsurv fit

Currently, only deviance residuals are supported.

## Usage

``` r
# S3 method for class 'ncvreg'
residuals(object, lambda, which = 1:length(object$lambda), drop = TRUE, ...)
```

## Arguments

- object:

  Object of class `ncvreg` or `ncvsurv`.

- lambda:

  Values of the regularization parameter at which residuals are
  requested (numeric vector). For values of lambda not in the sequence
  of fitted models, linear interpolation is used.

- which:

  Index of the penalty parameter at which residuals are requested
  (default = all indices). If `lambda` is specified, this take
  precedence over `which`.

- drop:

  By default, if a single value of lambda is supplied, a vector of
  residuals is returned (logical; default=`TRUE`). Set `drop=FALSE` if
  you wish to have the function always return a matrix (see
  [`drop()`](https://rdrr.io/r/base/drop.html)).

- ...:

  Not used.

## Examples

``` r
data(Prostate)
X <- Prostate$X
y <- Prostate$y
fit <- ncvreg(X, y)
residuals(fit)[1:5, 1:5]
#>     0.84343   0.78658   0.73357   0.68413   0.63802
#> 1 -2.909170 -2.768833 -2.637955 -2.515898 -2.402066
#> 2 -2.640906 -2.470432 -2.311447 -2.163178 -2.024901
#> 3 -2.640906 -2.505586 -2.379387 -2.261693 -2.151932
#> 4 -2.640906 -2.455181 -2.281973 -2.120439 -1.969792
#> 5 -2.106823 -2.063294 -2.022698 -1.984838 -1.949530
head(residuals(fit, lambda=0.1))
#>           1           2           3           4           5           6 
#> -1.18876659 -1.04838206 -0.90774951 -0.91188292 -1.52674317 -0.03250198 
```
