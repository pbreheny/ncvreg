[![version](http://www.r-pkg.org/badges/version/ncvreg)](https://cran.r-project.org/package=ncvreg)
[![downloads](http://cranlogs.r-pkg.org/badges/ncvreg)](https://cran.r-project.org/package=ncvreg)
[![codecov.io](https://codecov.io/github/pbreheny/ncvreg/coverage.svg?branch=master)](https://codecov.io/github/pbreheny/ncvreg?branch=master)
[![Travis build
status](https://travis-ci.org/pbreheny/breheny.svg?branch=master)](https://travis-ci.org/pbreheny/breheny)

# Regularization paths for MCP and SCAD penalized regression models

`ncvreg` is an R package for fitting regularization paths for linear
regression, GLM, and Cox regression models using lasso or nonconvex
penalties, in particular the minimax concave penalty (MCP) and smoothly
clipped absolute deviation (SCAD) penalty, with options for additional
L<sub>2</sub> penalties (the "elastic net" idea). Utilities for carrying
out cross-validation as well as post-fitting visualization,
summarization, inference, and prediction are also provided.

To learn more about the usage and syntax of `ncvreg`, see the [vignette](http://pbreheny.github.io/ncvreg/articles/getting-started.html) and follow the links under "Learn more".  To learn more about the algorithms used by `ncvreg`, see the original article:

  - [Breheny P and Huang J (2011). Coordinate descent algorithms for
    nonconvex penalized regression, with applications to biological
    feature selection. *Annals of Applied Statistics*, 5:
    232–253](http://myweb.uiowa.edu/pbreheny/pdf/Breheny2011.pdf)

For more about the marginal false discovery rate idea used for
post-selection inference, see

  - [Breheny P (2019). Marginal false discovery rates for penalized
    regression models. *Biostatistics*, **20**:
    299-314](https://dx.doi.org/10.1093/biostatistics/kxy004)

## Installation

`ncvreg` is on CRAN, so it can be installed via:

``` r
install.packages("ncvreg")
```
