# Plot coefficients from a ncvreg object

Produces a plot of the coefficient paths for a fitted `ncvreg` object.

## Usage

``` r
# S3 method for class 'ncvreg'
plot(x, alpha = 1, log.l = FALSE, shade = TRUE, col, ...)
```

## Arguments

- x:

  Fitted `"ncvreg"` model.

- alpha:

  Controls alpha-blending, helpful when the number of features is large.
  Default is alpha=1.

- log.l:

  Should horizontal axis be on the log scale? Default is FALSE.

- shade:

  Should nonconvex region be shaded? Default is TRUE.

- col:

  Vector of colors for coefficient lines. By default, evenly spaced
  colors are selected automatically.

- ...:

  Other graphical parameters to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

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
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y)
plot(fit)

plot(fit, col="black")

plot(fit, log=TRUE)

fit <- ncvreg(Prostate$X, Prostate$y, penalty.factor=rep(c(1, 1, 1, Inf), 2))
plot(fit, col=c('red', 'black', 'green'))  # Recycled among nonzero paths
```
