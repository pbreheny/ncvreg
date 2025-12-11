# ncvreg: Regularization Paths for SCAD and MCP Penalized Regression Models

Fits regularization paths for linear regression, GLM, and Cox regression
models using lasso or nonconvex penalties, in particular the minimax
concave penalty (MCP) and smoothly clipped absolute deviation (SCAD)
penalty, with options for additional L2 penalties (the "elastic net"
idea). Utilities for carrying out cross-validation as well as
post-fitting visualization, summarization, inference, and prediction are
also provided. For more information, see Breheny and Huang (2011)
[doi:10.1214/10-AOAS388](https://doi.org/10.1214/10-AOAS388) or visit
the ncvreg homepage <https://pbreheny.github.io/ncvreg/>.

## References

Breheny P and Huang J. (2011) Coordinate descent algorithms for
nonconvex penalized regression, with applications to biological feature
selection. *Annals of Applied Statistics*, **5**: 232-253.

## See also

Useful links:

- <https://pbreheny.github.io/ncvreg/>

- <https://github.com/pbreheny/ncvreg>

- Report bugs at <https://github.com/pbreheny/ncvreg/issues>

## Author

**Maintainer**: Patrick Breheny <patrick-breheny@uiowa.edu>
([ORCID](https://orcid.org/0000-0002-0650-1119))

Authors:

- Ryan Miller ([ORCID](https://orcid.org/0000-0003-0446-9992))

- Logan Harris ([ORCID](https://orcid.org/0000-0001-8562-9534))

## Examples

``` r
vignette("getting-started", package="ncvreg")
#> Warning: vignette ‘getting-started’ not found
```
