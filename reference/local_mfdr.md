# Estimate local mFDR for all features

`local_mfdr()` is called by
[`summary.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/summary.ncvreg.md),
which typically offers a more convenient interface to users. If,
however, you are working with local mfdrs programmatically rather than
interactively, you probably want to use `local_mfdr()`, which skips the
sorting, filtering, and print formatting of
[`summary.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/summary.ncvreg.md).

## Usage

``` r
local_mfdr(
  fit,
  lambda,
  X = NULL,
  y = NULL,
  method = c("ashr", "kernel"),
  sigma,
  ...
)
```

## Arguments

- fit:

  A fitted `ncvreg` or `ncvsurv` object.

- lambda:

  The value of lambda at which inference should be carried out.

- X, y:

  The design matrix and response used to fit the model; in most cases,
  it is not necessary to provide `X` and `y` as they are returned by
  `ncvreg`, but see the `returnX` argument in
  [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md).

- method:

  What method should be used to calculate the local fdr? Options are
  `ashr` (which tends to be more accurate) and `kernel` (which requires
  no additional packages). The default is to use `ashr` if the package
  is installed.

- sigma:

  For linear regression models, users can supply an estimate of the
  residual standard deviation. The default is to use RSS / DF, where
  degrees of freedom are approximated using the number of nonzero
  coefficients.

- ...:

  Additional arguments to
  [`ashr::ash()`](https://rdrr.io/pkg/ashr/man/ash.html) if using
  `method='ashr'`.

## Value

If all features are penalized, then the object returns a data frame with
one row per feature and four columns:

- `Estimate`: The coefficient estimate from the penalized regression fit

- `z`: A test statistic that approximately follows a standard normal
  distribution under the null hypothesis that the feature is marginally
  independent of the outcome

- `mfdr`: The estimated marginal local false discovery rate

- `Selected`: Features with nonzero coefficient estimates are given an
  asterisk

If some features are penalized and others are not, then a list is
returned with two elements: `pen.vars`, which consists of the data frame
described above, and `unpen.vars`, a data frame with four columns:
`Estimate`, `SE`, `Statistic`, and `p.value`. The standard errors and
p-values are based on a classical `lm`/`glm`/`coxph` model using the
effect of the penalized features as an offset.

## See also

[`summary.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/summary.ncvreg.md)

## Examples

``` r
# Linear regression
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y)
local_mfdr(fit, 0.1)
#>          Estimate          z         mfdr Selected
#> lcavol  0.5387509  8.8807209 4.956187e-16        *
#> lweight 0.6382717  3.9554154 1.172927e-02        *
#> age     0.0000000 -1.0067558 9.454421e-01         
#> lbph    0.0000000  1.2720071 9.275879e-01         
#> svi     0.6102800  3.7616622 2.430660e-02        *
#> lcp     0.0000000 -0.1448424 9.660606e-01         
#> gleason 0.0000000  0.7670990 9.554227e-01         
#> pgg45   0.0000000  0.9391426 9.487372e-01         

fit <- ncvreg(Prostate$X, Prostate$y, penalty.factor=rep(0:1, each=4))
local_mfdr(fit, 0.1)
#> $pen.vars
#>          Estimate          z       mfdr Selected
#> svi     0.6965651  4.1620116 0.00927594        *
#> lcp     0.0000000 -0.3000954 0.97982329         
#> gleason 0.0000000  0.8748507 0.97194659         
#> pgg45   0.0000000  1.0057234 0.96838914         
#> 
#> $unpen.vars
#>            Estimate  std.error statistic      p.value
#> lcavol   0.54734893 0.06407513  8.542219 2.654740e-13
#> lweight  0.58932799 0.19646341  2.999457 3.479201e-03
#> age     -0.01641325 0.01061286 -1.546875 1.253265e-01
#> lbph     0.10049326 0.05671452  1.772200 7.967190e-02
#> 

# Logistic regression
data(Heart)
X <- Heart$X
y <- Heart$y
fit <- ncvreg(X, y, family='binomial')
local_mfdr(fit, 0.1)
#>             Estimate         z         mfdr Selected
#> sbp       0.00000000 2.1766234 4.935585e-01         
#> tobacco   0.00000000 4.1575377 1.883630e-03         
#> ldl       0.00000000 4.2061346 1.542018e-03         
#> adiposity 0.00000000 2.5179307 3.043436e-01         
#> famhist   0.05505671 4.7878991 1.190670e-04        *
#> typea     0.00000000 2.8422109 1.550916e-01         
#> obesity   0.00000000 0.7961026 8.835105e-01         
#> alcohol   0.00000000 0.9210390 8.720085e-01         
#> age       0.03577763 8.1799158 1.199132e-13        *

# Cox regression
data(Lung)
X <- Lung$X
y <- Lung$y
fit <- ncvsurv(X, y)
local_mfdr(fit, 0.1)
#>             Estimate          z         mfdr Selected
#> trt       0.03737689  1.3316047 9.462748e-01        *
#> karno    -0.03196343 -6.3518252 1.081598e-07        *
#> diagtime  0.00000000  0.2972526 9.761319e-01         
#> age       0.00000000 -0.5655773 9.732810e-01         
#> prior     0.00000000  0.2907796 9.761762e-01         
#> squamous -0.64313783 -3.1199729 2.476091e-01        *
#> small     0.00000000  0.2041520 9.766697e-01         
#> adeno     0.28016650  2.0507762 8.392153e-01        *
#> large    -0.20731421 -1.7987677 8.944907e-01        *
```
