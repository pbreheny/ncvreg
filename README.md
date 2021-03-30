[![GitHub version](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/pbreheny/ncvreg/master/.version.json&style=flat&logo=github)](https://github.com/pbreheny/ncvreg)
[![CRAN version](https://img.shields.io/cran/v/ncvreg?logo=R)](https://cran.r-project.org/package=ncvreg)
[![downloads](https://cranlogs.r-pkg.org/badges/ncvreg)](https://cran.r-project.org/package=ncvreg)
[![Travis build status](https://travis-ci.org/pbreheny/ncvreg.svg?branch=master)](https://travis-ci.org/pbreheny/ncvreg)
[![codecov.io](https://codecov.io/github/pbreheny/ncvreg/coverage.svg?branch=master)](https://codecov.io/github/pbreheny/ncvreg?branch=master)

# [Regularization paths for MCP and SCAD penalized regression models](https://pbreheny.github.io/ncvreg/)

`ncvreg` is an R package for fitting regularization paths for linear
regression, GLM, and Cox regression models using lasso or nonconvex
penalties, in particular the minimax concave penalty (MCP) and smoothly
clipped absolute deviation (SCAD) penalty, with options for additional
L<sub>2</sub> penalties (the "elastic net" idea). Utilities for carrying
out cross-validation as well as post-fitting visualization,
summarization, inference, and prediction are also provided.

* To get started using `ncvreg`, see the ["getting started" vignette](https://pbreheny.github.io/ncvreg/articles/getting-started.html)
* To learn more, follow the links under "Learn more" at the [ncvreg website](https://pbreheny.github.io/ncvreg/)
* For details on the algorithms used by `ncvreg`, see the original article: [Breheny P and Huang J (2011). Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection. *Annals of Applied Statistics*, 5: 232–253](https://myweb.uiowa.edu/pbreheny/pdf/Breheny2011.pdf)
* For more about the marginal false discovery rate idea used for
post-selection inference, see [Breheny P (2019). Marginal false discovery rates for penalized regression models. *Biostatistics*, **20**: 299-314](https://dx.doi.org/10.1093/biostatistics/kxy004)

## Installation

To install the latest release version from CRAN:

```r
install.packages("ncvreg")
```

To install the latest development version from GitHub:

```r
remotes::install_github("pbreheny/ncvreg")
```

