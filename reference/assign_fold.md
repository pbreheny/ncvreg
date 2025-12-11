# Assign folds for cross-validation

If y has only two unique values, fold assignments are chosen so that the
balance between outcomes is the same in each fold. This is useful for
logistic regression and time-to-event data (to balance the fraction of
observations that are censored).

## Usage

``` r
assign_fold(y, folds, seed)
```

## Arguments

- y:

  Either (i) the vector of outcomes or (ii) a vector such as `1:n`

- folds:

  Number of folds

- seed:

  (optional) set a seed for reproducibility

## Value

A vector of integers indicating fold assignments

## See also

`[cv.ncvreg()]`

## Examples

``` r
assign_fold(rnorm(11), 2)
#>  [1] 1 1 2 1 1 2 2 2 1 2 1
assign_fold(1:41, 7)
#>  [1] 3 1 3 3 2 4 4 2 5 1 5 1 5 3 2 2 5 6 4 6 1 7 6 7 5 3 4 7 1 2 3 7 1 4 7 5 6 6
#> [39] 2 4 6
assign_fold(1:41, 7) |> table()
#> 
#> 1 2 3 4 5 6 7 
#> 6 6 6 6 6 6 5 
data(Heart)
assign_fold(Heart$y, 7) |> head()
#> [1] 2 5 1 4 2 3
assign_fold(Heart$y, 7) |> table(Heart$y)
#>    
#>      0  1
#>   1 43 23
#>   2 43 23
#>   3 43 23
#>   4 43 23
#>   5 43 23
#>   6 43 23
#>   7 44 22
```
