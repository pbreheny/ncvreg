# Permute residuals for a fitted ncvreg model

Fits multiple penalized regression models in which the residuals are
randomly permuted, thereby allowing estimation of the marginal false
discovery rate.

## Usage

``` r
permres(fit, ...)

# S3 method for class 'ncvreg'
permres(fit, lambda, N = 10, seed, trace = FALSE, ...)
```

## Arguments

- fit:

  A fitted ncvreg model, as produced by
  [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md).
  To use with `permres`, the model must be fit using the `returnX=TRUE`
  option.

- ...:

  Not used.

- lambda:

  The regularization parameter to use for estimating residuals. Unlike
  [`perm.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/perm.ncvreg.md),
  `permres()` calculates EF and mFDR for a specific `lambda` value, not
  an entire path. As a result, it runs much faster.

- N:

  The number of permutation replications. Default is 10.

- seed:

  You may set the seed of the random number generator in order to obtain
  reproducible results.

- trace:

  If set to TRUE, perm.ncvreg will inform the user of its progress by
  announcing the beginning of each permutation fit. Default is FALSE.

## Value

A list with the following components:

- EF:

  The number of variables selected at each value of `lambda`, averaged
  over the permutation fits.

- S:

  The actual number of selected variables for the non-permuted data.

- mFDR:

  The estimated marginal false discovery rate (`EF/S`).

- loss:

  The loss/deviance, averaged over the permutation fits. This is an
  estimate of the explanatory power of the model under null conditions,
  and can be used to adjust the loss of the fitted model in a manner
  akin to the idea of an adjusted R-squared in classical regression.

## Details

The function fits a penalized regression model to the actual data, then
repeats the process `N` times with a permuted version of the response
vector. This allows estimation of the expected number of variables
included by chance for each value of `lambda`. The ratio of this
expected quantity to the number of selected variables using the actual
(non-permuted) response is called the marginal false discovery rate
(mFDR).

## See also

[`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md),
\`[`mfdr()`](https://pbreheny.github.io/ncvreg/reference/mfdr.md),
[`perm.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/perm.ncvreg.md)

## Author

Patrick Breheny <patrick-breheny@uiowa.edu>

## Examples

``` r
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y, N=50)
permres(fit, lambda=0.15)
#> $EF
#> [1] 0.2
#> 
#> $S
#> 0.1500 
#>      3 
#> 
#> $mFDR
#>     0.1500 
#> 0.06666667 
#> 
#> $loss
#> [1] 49.14291
#> 
```
