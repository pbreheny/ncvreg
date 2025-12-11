# Projection base test statistics and intervals

Constructs projection based test statistics that can be used to control
FDR along with intervals for a penalized regression model.

## Usage

``` r
intervals(
  fit,
  lambda,
  sigma,
  level = 0.95,
  posterior = TRUE,
  relaxed = FALSE,
  adjust_projection = FALSE,
  X = NULL
)
```

## Arguments

- fit:

  An optional fit of class `ncvreg` or `cv.ncvreg`. If supplied, `X`
  should only be supplied if `fit` does not contain it.

- lambda:

  The penalty at which the tests and intervals are to be constructed. If
  left unspecified, will be selected using cross validation.

- sigma:

  Standard deviation estimate used to compute test statistic and
  intervals. If left unspecified (default) it will be estimated using
  the recommendation from Reid et al. (2016)

- level:

  the confidence level required.

- posterior:

  whether the intervals returned should be posterior intervals (default)
  or debiased intervals (if `FALSE`). Posterior intervals are
  constructed from distributions where the coefficient estimates are the
  is the posterior mode. Debiased intervals are constructed around the
  estimates.

- relaxed:

  whether the relaxed lasso based statistic / intervals should be used.
  Default is `FALSE` in which case PIPE based intervals are constructed
  (recommended). This affects the estimate.

- adjust_projection:

  whether a Local Quadratic Approximation should be used in the
  projection for determining variance of the PIPE based test statistic.
  Default is `FALSE` as more research has been done without this
  adjustment, however without this adjustment, statistics and intervals
  may be over conservative in the presence of correlation. This affects
  the SE.

- X:

  The original design matrix supplied to `fit`. Required if `fit` does
  not contain `X`.

## Value

An `data.frame` containing the following columns:

- variable:

  `colnames(X)`

- coef:

  The original estimates at the specified parameters (`lambda`, `gamma`,
  `alpha`)

- estimate:

  The debiased estimates.

- SE:

  The standard errors.

- t:

  The PIPE / Relaxed Lasso / LQA test statistics

- lower:

  Interval lower bounds

- upper:

  Intervals upper bounds

- p.value:

  The unadjusted p-value

- p.adjust:

  The Benhamini and Hochberg corrected p-value

- penalty:

  The penalty used.

- lambda:

  The lambda value the test statistics and intervals were constructed
  at.

- gamma:

  The gamma value the test statistics and intervals were constructed at
  (for MCP/SCAD).

- alpha:

  The alpha value the test statistics and intervals were constructed at.

- level:

  The confidence level set for interval construction.

- sigma:

  The standard deviation used for constructing the test statistis and
  intervals.

## Details

The function constructs test statistics and intervals based off an
approximate projection onto the column space of the active features. The
test statistic can be used to control FDR and the intervals generally
have good coverage. However, both tend to be conservative with the
introduction of correlation.

The intervals produced can either be biased (like the point estimates)
or debiased by setting the parameter `posterior` accordingly. The
resulting behavior is quite different. See references for more details.

## References

Harris L and Breheny P. (2025) A new perspective on high dimensional
confidence intervals. *arXiv preprint*, arXiv:2508.03504.
<https://arxiv.org/abs/2508.03504>

Harris L and Breheny P. (2025) Alternative Likelihood Approximations for
High-Dimensional Intervals for Lasso. *arXiv preprint*,
arXiv:2509.14971. <https://arxiv.org/abs/2509.14971>

Dai B. (2019) Projection-based inference and model selection for
penalized regression. PhD dissertation, University of Iowa, Iowa City,
IA. [doi:10.17077/etd.005250](https://doi.org/10.17077/etd.005250)

## Author

Logan Harris, Patrick Breheny, and Biyue Dai

## Examples

``` r
# Linear regression (SCAD-Net penalty, PIPE intervals, pass ncvreg object)
fit <- ncvreg(Prostate$X, Prostate$y, penalty = "SCAD", alpha = 0.9)
intervals(fit) |> head()
#>         variable        coef    estimate         SE         t       lower
#> lcavol    lcavol  0.56705170  0.56823231 0.08536950  6.656151  0.39993486
#> lweight  lweight  0.61396251  0.61523624 0.19734546  3.117560  0.22673738
#> age          age -0.02075306 -0.02079495 0.01091667 -1.904881 -0.04210383
#> lbph        lbph  0.09698225  0.09717897 0.05726386  1.697039 -0.01541831
#> svi          svi  0.74938522  0.75104331 0.23685484  3.170901  0.28524770
#> lcp          lcp -0.10234574 -0.10256268 0.08884987 -1.154337 -0.27616791
#>                upper      p.value     p.adjust penalty    lambda gamma alpha
#> lcavol  0.7342398444 2.810907e-11 2.248726e-10    SCAD 0.0201901   3.7   0.9
#> lweight 1.0003747644 1.823550e-03 4.862801e-03    SCAD 0.0201901   3.7   0.9
#> age     0.0008049362 5.679554e-02 1.135911e-01    SCAD 0.0201901   3.7   0.9
#> lbph    0.2089392044 8.968937e-02 1.435030e-01    SCAD 0.0201901   3.7   0.9
#> svi     1.2132756360 1.519668e-03 4.862801e-03    SCAD 0.0201901   3.7   0.9
#> lcp     0.0714558667 2.483621e-01 2.838424e-01    SCAD 0.0201901   3.7   0.9
#>         level    sigma
#> lcavol   0.95 0.692084
#> lweight  0.95 0.692084
#> age      0.95 0.692084
#> lbph     0.95 0.692084
#> svi      0.95 0.692084
#> lcp      0.95 0.692084

# Logistic regression (lasso penalty, LQA intervals, pass cv.ncvreg object) 
data(Heart)
cv_fit <- cv.ncvreg(Heart$X, Heart$y, family="binomial", penalty = "lasso")
intervals(cv_fit, adjust_projection = TRUE) |> head()
#>            variable        coef    estimate         SE         t        lower
#> sbp             sbp 0.004591841 0.006691921 0.00542863 1.2327090 -0.004042523
#> tobacco     tobacco 0.072039912 0.081846321 0.02488849 3.2885212  0.023402587
#> ldl             ldl 0.153354310 0.175073685 0.05424203 3.2276386  0.047449143
#> adiposity adiposity 0.000000000 0.004350699 0.01881224 0.2312696 -0.029110113
#> famhist     famhist 0.828918109 0.917515358 0.22125877 4.1467977  0.395299160
#> typea         typea 0.031161890 0.036060712 0.01169760 3.0827432  0.008391284
#>                upper      p.value     p.adjust penalty      lambda gamma alpha
#> sbp       0.01537976 2.176844e-01 0.3244134078   lasso 0.008236943     3     1
#> tobacco   0.12081978 1.007152e-03 0.0028083724   lasso 0.008236943     3     1
#> ldl       0.25967645 1.248165e-03 0.0028083724   lasso 0.008236943     3     1
#> adiposity 0.03657256 8.171053e-01 0.8810252658   lasso 0.008236943     3     1
#> famhist   1.26256560 3.371576e-05 0.0001517209   lasso 0.008236943     3     1
#> typea     0.05409455 2.051020e-03 0.0036918367   lasso 0.008236943     3     1
#>           level sigma
#> sbp        0.95    NA
#> tobacco    0.95    NA
#> ldl        0.95    NA
#> adiposity  0.95    NA
#> famhist    0.95    NA
#> typea      0.95    NA
```
