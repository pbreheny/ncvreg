[![version](http://www.r-pkg.org/badges/version/ncvreg)](https://cran.r-project.org/package=ncvreg)
[![downloads](http://cranlogs.r-pkg.org/badges/ncvreg)](https://cran.r-project.org/package=ncvreg)
[![codecov.io](https://codecov.io/github/pbreheny/ncvreg/coverage.svg?branch=master)](https://codecov.io/github/pbreheny/ncvreg?branch=master)

# Regularization Paths for SCAD and MCP Penalized Regression Models

`ncvreg` fits regularization paths for linear regression, GLM, and Cox regression models using lasso or nonconvex penalties, in particular the minimax concave penalty (MCP) and smoothly clipped absolute deviation (SCAD) penalty, with options for additional L<sub>2</sub> penalties (the "elastic net" idea).  Utilities for carrying out cross-validation as well as post-fitting visualization, summarization, inference, and prediction are also provided.

## Basic Usage

The basic usage of ncvreg is as follows:

```r
fit <- ncvreg(X, y)
```

The default penalty here is the minimax concave penalty (MCP), but SCAD and lasso penalties are also available.  This produces a path of coefficients, which we can plot with

```r
plot(fit)
```

<p align="center">
<img alt="img" width=480 src="http://pbreheny.github.io/ncvreg/index_files/figure-html/plot-1.png">
</p>

Notice that variables enter the model one at a time, and that at any given value of `lambda`, several coefficients are zero.  The `summary` method can be used for post-selection summarization and inference:

```{r summary}
summary(fit, lambda=0.05)

# MCP-penalized linear regression with n=97, p=8
# At lambda=0.0500:
# -------------------------------------------------
#   Nonzero coefficients: 6
#   Expected nonzero coefficients: 2.51
#   Average mfdr (6 features)    : 0.418
```

`summary(fit)` also returns the following table:

|        |   Estimate|         z|      mfdr|
|:-------|----------:|---------:|---------:|
|lcavol  |  0.5317899|  8.880429| 0.0000000|
|svi     |  0.6725610|  3.945052| 0.0018967|
|lweight |  0.6038969|  3.665874| 0.0050683|
|lbph    |  0.0887456|  1.928241| 0.4998035|
|age     | -0.0153092| -1.788334| 1.0000000|
|pgg45   |  0.0016804|  1.159772| 1.0000000|

In this case, it would appear that `lcavol`, `svi`, and `lweight` are clearly associated with the response, even after adjusting for the other variables in the model, while `lbph`, `age`, and `pgg45` may be false positives included simply by chance.

Typically, one would carry out cross-validation for the purposes of assessing the predictive accuracy of the model at various values of `lambda`:

```{r cvplot, h=4, w=6, cache=TRUE}
cvfit <- cv.ncvreg(X, y)
plot(cvfit)
```

<p align="center">
<img alt="img" width=480 src="http://pbreheny.github.io/ncvreg/index_files/figure-html/cvplot-1.png">
</p>

At this point, `coef(cvfit)` will return the coefficients at the value of `lambda` minimizing the cross-validation error.  Likewise,

```r
predict(cvfit, X=head(X))
```

will return predictions for that model, while

```r
predict(cvfit, type="nvars")
```

will return the number of nonzero coefficients.  Note that the original fit (to the full data set) is returned as `cvfit$fit`; it is not necessary to call both `ncvreg` and `cv.ncvreg` to analyze a data set.  For example, `plot(cvfit$fit)` will produce the same coefficient path plot as `plot(fit)` above.

## Documentation and Citation

For more on the usage and syntax of `ncvreg`, see the [ncvreg homepage](http://pbreheny.github.io/ncvreg).

For more on the algorithms used by `ncvreg`, see the original article:

* [Breheny P and Huang J (2011).  Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection.  *Annals of Applied Statistics*, 5: 232--253](http://myweb.uiowa.edu/pbreheny/pdf/Breheny2011.pdf)

For more about the marginal false discovery rate idea used for post-selection inference, see

* [Breheny P (to appear).  Marginal false discovery rates for penalized regression models.  *Biostatistics*](https://arxiv.org/pdf/1607.05636)

## Installation

* To install the latest release version from CRAN: `install.packages("ncvreg")`
* To install the latest development version from GitHub: `devtools::install_github("pbreheny/ncvreg")`
