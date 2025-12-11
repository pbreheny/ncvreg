# Summary method for ncvreg objects

Inferential summaries for `ncvreg` and `ncvsurv` objects based on local
marginal false discovery rates.

## Usage

``` r
# S3 method for class 'ncvreg'
summary(object, lambda, which, number, cutoff, sort = TRUE, sigma, ...)

# S3 method for class 'summary.ncvreg'
print(x, digits, ...)
```

## Arguments

- object:

  An `ncvreg` or `ncvsurv` object.

- lambda:

  The regularization parameter value at which inference should be
  reported.

- which:

  Alternatively, `lambda` may be specified by index; `which=10` means:
  report inference for the 10th value of `lambda` along the
  regularization path. If both `lambda` and `which` are specified,
  `lambda` takes precedence.

- number:

  By default, `summary` will provide an inferential summary for each
  variable that has been selected (i.e., each variable with a nonzero
  coefficient). Specifying `number=5`, for example, means that the
  summary table will include the 5 features with the lowest mfdr values,
  regardless of whether they were selected. To see all features,
  `number=Inf`.

- cutoff:

  Alternatively, specifying for example `cutoff=0.3` will report
  inference for all features with mfdr under 30%. If both `number` and
  `cutoff` are specified, the intersection between both sets of features
  is reported.

- sort:

  Should the results be sorted by `mfdr`? (default: TRUE)

- sigma:

  For linear regression models, users can supply an estimate of the
  residual standard deviation. The default is to use RSS / DF, where
  degrees of freedom are approximated using the number of nonzero
  coefficients.

- ...:

  Further arguments; in particular, if you have set `returnX=FALSE`, you
  will need to supply `X` and `y` in order to calculate local mFDRs.

- x:

  A `summary.ncvreg` object.

- digits:

  Number of digits past the decimal point to print out. Can be a vector
  specifying different display digits for each of the five non-integer
  printed values.

## Value

An object with S3 class `summary.ncvreg`. The class has its own print
method and contains the following list elements:

- penalty:

  The penalty used by `ncvreg` or `ncvsurv`

- model:

  Either `"linear"`, `"logistic"`, or `"Cox"`.

- n:

  Number of instances.

- p:

  Number of regression coefficients (not including the intercept).

- lambda:

  The `lambda` value at which inference is being reported.

- nvars:

  The number of nonzero coefficients (again, not including the
  intercept) at that value of `lambda`.

- table:

  A table containing estimates, normalized test statistics (z), and an
  estimate of the local mfdr for each coefficient. The mfdr may be
  loosely interpreted, in an empirical Bayes sense, as the probability
  that the given feature is null.

- unpen.table:

  If there are any unpenalized coefficients, a separate inferential
  summary is given for them. Currently, this is based on
  `lm`/`glm`/`coxph` using the penalized coefficients to provide an
  offset. This is useful and more or less accurate, but not ideal; we
  hope to improve the inferential methods for unpenalized variables in
  the future.

## See also

[`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md),
[`cv.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.md),
[`plot.cv.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/plot.cv.ncvreg.md),
[`local_mfdr()`](https://pbreheny.github.io/ncvreg/reference/local_mfdr.md)

## Author

Patrick Breheny <patrick-breheny@uiowa.edu>

## Examples

``` r
# Linear regression --------------------------------------------------
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y)
summary(fit, lambda=0.08)
#> MCP-penalized linear regression with n=97, p=8
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Nonzero coefficients         :   5
#>   Expected nonzero coefficients:   1.90
#>   Average mfdr (5 features)    :   0.381
#> 
#>           Estimate      z     mfdr Selected
#> lcavol   5.263e-01  8.618  < 1e-04        *
#> svi      6.724e-01  3.867 0.016288        *
#> lweight  6.410e-01  3.815 0.019837        *
#> lbph     1.322e-02  1.295 0.926736        *
#> age     -6.203e-05 -1.095 0.941368        *

# Logistic regression ------------------------------------------------
data(Heart)
fit <- ncvreg(Heart$X, Heart$y, family="binomial")
summary(fit, lambda=0.05)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0500:
#> -------------------------------------------------
#>   Nonzero coefficients         :   5
#>   Expected nonzero coefficients:   0.18
#>   Average mfdr (5 features)    :   0.036
#> 
#>         Estimate     z      mfdr Selected
#> age     0.047479 6.645   < 1e-04        *
#> famhist 0.611522 4.332 0.0010983        *
#> tobacco 0.037272 3.465 0.0308975        *
#> ldl     0.075045 3.385 0.0403379        *
#> typea   0.009799 3.057 0.1077203        *

# Cox regression -----------------------------------------------------
data(Lung)
fit <- ncvsurv(Lung$X, Lung$y)
summary(fit, lambda=0.1)
#> MCP-penalized Cox regression with n=137, p=8
#> At lambda=0.1000:
#> -------------------------------------------------
#>   Nonzero coefficients         :   5
#>   Expected nonzero coefficients:   2.93
#>   Average mfdr (5 features)    :   0.586
#> 
#>          Estimate      z    mfdr Selected
#> karno    -0.03196 -6.352 < 1e-04        *
#> squamous -0.64314 -3.120 0.24761        *
#> adeno     0.28017  2.051 0.83922        *
#> large    -0.20731 -1.799 0.89449        *
#> trt       0.03738  1.332 0.94627        *

# Options ------------------------------------------------------------
fit <- ncvreg(Heart$X, Heart$y, family="binomial")
summary(fit, lambda=0.08, number=3)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Features satisfying criteria       : 3
#>   Average mfdr among chosen features : 0.0026
#> 
#>         Estimate     z       mfdr Selected
#> age     0.041921 7.771    < 1e-04        *
#> famhist 0.277822 4.627 0.00029409        *
#> ldl     0.009273 3.847 0.00751773        *
summary(fit, lambda=0.08, number=Inf)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Features satisfying criteria       : 9
#>   Average mfdr among chosen features : 0.377
#> 
#>           Estimate      z       mfdr Selected
#> age       0.041921 7.7712    < 1e-04        *
#> famhist   0.277822 4.6274 0.00029409        *
#> ldl       0.009273 3.8473 0.00751773        *
#> tobacco   0.005783 3.7881 0.00938405        *
#> typea     0.000000 2.9383 0.13986436         
#> adiposity 0.000000 1.8121 0.70188341         
#> sbp       0.000000 1.7780 0.71453183         
#> alcohol   0.000000 0.7807 0.89965322         
#> obesity   0.000000 0.4339 0.91712823         
summary(fit, lambda=0.08, cutoff=0.5)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Features satisfying criteria       : 5
#>   Average mfdr among chosen features : 0.0314
#> 
#>         Estimate     z       mfdr Selected
#> age     0.041921 7.771    < 1e-04        *
#> famhist 0.277822 4.627 0.00029409        *
#> ldl     0.009273 3.847 0.00751773        *
#> tobacco 0.005783 3.788 0.00938405        *
#> typea   0.000000 2.938 0.13986436         
summary(fit, lambda=0.08, number=3, cutoff=0.5)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Features satisfying criteria       : 3
#>   Average mfdr among chosen features : 0.0026
#> 
#>         Estimate     z       mfdr Selected
#> age     0.041921 7.771    < 1e-04        *
#> famhist 0.277822 4.627 0.00029409        *
#> ldl     0.009273 3.847 0.00751773        *
summary(fit, lambda=0.08, number=5, cutoff=0.1)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Features satisfying criteria       : 4
#>   Average mfdr among chosen features : 0.0043
#> 
#>         Estimate     z       mfdr Selected
#> age     0.041921 7.771    < 1e-04        *
#> famhist 0.277822 4.627 0.00029409        *
#> ldl     0.009273 3.847 0.00751773        *
#> tobacco 0.005783 3.788 0.00938405        *
summary(fit, lambda=0.08, number=Inf, sort=FALSE)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Features satisfying criteria       : 9
#>   Average mfdr among chosen features : 0.377
#> 
#>           Estimate      z       mfdr Selected
#> sbp       0.000000 1.7780 0.71453183         
#> tobacco   0.005783 3.7881 0.00938405        *
#> ldl       0.009273 3.8473 0.00751773        *
#> adiposity 0.000000 1.8121 0.70188341         
#> famhist   0.277822 4.6274 0.00029409        *
#> typea     0.000000 2.9383 0.13986436         
#> obesity   0.000000 0.4339 0.91712823         
#> alcohol   0.000000 0.7807 0.89965322         
#> age       0.041921 7.7712    < 1e-04        *
summary(fit, lambda=0.08, number=3, cutoff=0.5, sort=FALSE)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Features satisfying criteria       : 3
#>   Average mfdr among chosen features : 0.0026
#> 
#>         Estimate     z       mfdr Selected
#> ldl     0.009273 3.847 0.00751773        *
#> famhist 0.277822 4.627 0.00029409        *
#> age     0.041921 7.771    < 1e-04        *

# If X and y are not returned with the fit, they must be supplied
fit <- ncvreg(Heart$X, Heart$y, family="binomial", returnX=FALSE)
summary(fit, X=Heart$X, y=Heart$y, lambda=0.08)
#> MCP-penalized logistic regression with n=462, p=9
#> At lambda=0.0800:
#> -------------------------------------------------
#>   Nonzero coefficients         :   4
#>   Expected nonzero coefficients:   0.02
#>   Average mfdr (4 features)    :   0.004
#> 
#>         Estimate     z       mfdr Selected
#> age     0.041921 7.771    < 1e-04        *
#> famhist 0.277822 4.627 0.00029409        *
#> ldl     0.009273 3.847 0.00751773        *
#> tobacco 0.005783 3.788 0.00938405        *
```
