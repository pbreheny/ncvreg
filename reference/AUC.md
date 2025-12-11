# AUC for cv.ncvsurv objects

Calculates the cross-validated AUC (concordance) from a `cv.ncvsurv`
object.

## Usage

``` r
# S3 method for class 'cv.ncvsurv'
AUC(obj, ...)
```

## Arguments

- obj:

  A `cv.ncvsurv` object. You must run
  [`cv.ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.md)
  with the option `returnY=TRUE` in order for `AUC()` to work.

- ...:

  For S3 method compatibility; not used

## Details

The area under the curve (AUC), or equivalently, the concordance
statistic (C), is calculated according to the procedure described in van
Houwelingen and Putter (2011). The function calls
[`survival::concordancefit()`](https://rdrr.io/pkg/survival/man/concordancefit.html),
except cross-validated linear predictors are used to guard against
overfitting. Thus, the values returned by `AUC.cv.ncvsurv()` will be
lower than those you would obtain with `concordancefit()` if you fit the
full (unpenalized) model.

## References

van Houwelingen H, Putter H (2011). Dynamic Prediction in Clinical
Survival Analysis. CRC Press.

## See also

[`cv.ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.md),
[`survival::concordancefit()`](https://rdrr.io/pkg/survival/man/concordancefit.html)

## Author

Patrick Breheny, Brandon Butcher, and Lawrence Hunsicker

## Examples

``` r
data(Lung)
X <- Lung$X
y <- Lung$y

cvfit <- cv.ncvsurv(X, y, returnY=TRUE)
head(AUC(cvfit))
#> [1] 0.5135734 0.6059746 0.6579964 0.6710586 0.6741254 0.6745797
lam <- cvfit$lambda
plot(lam, AUC(cvfit), xlim=rev(range(lam)), lwd=3, type='l',
     las=1, xlab=expression(lambda), ylab='AUC')
```
