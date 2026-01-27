# Plot survival curve for ncvsurv model

Plot survival curve for a model that has been fit using
[`ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/ncvsurv.md)
followed by a prediction of the survival function using
[`predict.ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/predict.ncvsurv.md).

## Usage

``` r
# S3 method for class 'ncvsurv.func'
plot(x, alpha = 1, ...)
```

## Arguments

- x:

  A `ncvsurv.func` object, which is returned by
  [`predict.ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/predict.ncvsurv.md)
  if `type='survival'` is specified. See examples.

- alpha:

  Controls alpha-blending (i.e., transparency). Useful if many
  overlapping lines are present.

- ...:

  Other graphical parameters to pass to `plot`

## See also

[`ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/ncvsurv.md),
[`predict.ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/predict.ncvsurv.md)

## Author

Patrick Breheny

## Examples

``` r
data(Lung)
X <- Lung$X
y <- Lung$y

fit <- ncvsurv(X, y)

# A single survival curve
S <- predict(fit, X[1,], type = "survival", lambda = 0.15)
plot(S, xlim = c(0, 200))


# Lots of survival curves
S <- predict(fit, X, type = "survival", lambda = 0.08)
plot(S, xlim = c(0, 200), alpha = 0.3)
```
