# Standardizes a design matrix

Accepts a design matrix and returns a standardized version of that
matrix (i.e., each column will have mean 0 and mean sum of squares equal
to 1).

## Usage

``` r
std(X, Xnew)
```

## Arguments

- X:

  A matrix (or object that can be coerced to a matrix, such as a data
  frame or numeric vector).

- Xnew:

  Optional. If supplied, `X` must be the output of `std()` and `Xnew` is
  to be standardized in the same way. See examples for why this might be
  useful.

## Value

The standardized design matrix, with the following attribues:

- center, scale:

  mean and standard deviation used to scale the columns

- nonsingular:

  A vector indicating which columns of the original design matrix were
  able to be standardized (constant columns cannot be standardized to
  have a standard deviation of 1)

## Details

This function centers and scales each column of `X` so that
\$\$\sum\_{i=1}^n x\_{ij}=0\$\$ and \$\$n^{-1} \sum\_{i=1}^n x\_{ij}^2 =
1\$\$ for all j. This is usually not necessary to call directly, as
**ncvreg** internally standardizes the design matrix, but inspection of
the standardized design matrix can sometimes be useful. This differs
from the base R function [`scale()`](https://rdrr.io/r/base/scale.html)
in two ways:

1.  [`scale()`](https://rdrr.io/r/base/scale.html) uses the sample
    standard deviation `sqrt(sum(x^2)/(n-1))`, while `std()` uses the
    root-mean-square standard deviation `sqrt(mean(sum(x^2))` without
    the \\n/(n-1)\\ correction

2.  `std` is faster.

## Examples

``` r
data(Prostate)
S <- std(Prostate$X)
apply(S, 2, sum)
#>        lcavol       lweight           age          lbph           svi 
#>  4.211909e-15  6.637225e-14  4.140785e-14  1.110223e-16 -2.886580e-15 
#>           lcp       gleason         pgg45 
#>  1.576864e-14 -4.218847e-15  4.982126e-15 
apply(S, 2, function(x) mean(x^2))
#>  lcavol lweight     age    lbph     svi     lcp gleason   pgg45 
#>       1       1       1       1       1       1       1       1 

# Standardizing new observations
X1 <- Prostate$X[1:90,]
X2 <- Prostate$X[91:97,]
S <- std(X1)
head(std(S, X2))
#>      lcavol    lweight        age       lbph        svi        lcp    gleason
#> 91 1.845841  1.1196037  0.5439694 -1.0472200 -0.4472136 -0.8399653 -1.0076612
#> 92 1.197790  0.1447897 -0.4241118  0.8380061  2.2360680 -0.8399653  0.3459136
#> 93 1.467844  0.6016472  0.5439694 -1.0472200  2.2360680  1.2129092  0.3459136
#> 94 2.367589  0.6487805 -2.7751661 -1.0472200  2.2360680  1.8552150  0.3459136
#> 95 1.537935 -0.5017477 -1.6687876 -1.0472200  2.2360680  2.0786918  0.3459136
#> 96 1.515337  0.3661621  0.5439694  0.9828411  2.2360680  1.3921070  0.3459136
#>         pgg45
#> 91 -0.8476607
#> 92 -0.3129215
#> 93  1.2912962
#> 94  0.5783106
#> 95 -0.4911679
#> 96  2.0042819
# Useful if you fit to a standardized X, but then get new obs:
y <- Prostate$y[1:90]
fit <- ncvreg(S, y)
predict(fit, std(S, X2), lambda=0.1)
#>       91       92       93       94       95       96       97 
#> 3.514077 2.938705 3.408371 3.813543 2.900628 3.445538 3.617198 
# Same as
predict(ncvreg(X1, y), X2, lambda=0.1)
#>       91       92       93       94       95       96       97 
#> 3.514077 2.938705 3.408371 3.813543 2.900628 3.445538 3.617198 
```
