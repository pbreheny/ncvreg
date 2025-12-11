# Risk factors associated with heart disease

Data from a subset of the Coronary Risk-Factor Study baseline survey,
carried out in rural South Africa.

## Usage

``` r
Heart
```

## Format

A list of two objects: `y` and `X`

- y:

  Coronary heart disease at baseline; 1=Yes 0=No

- X:

  A matrix with 462 observations (rows) and 9 predictor variables
  (columns). The remainder of this list describes the columns of `X`

- sbp:

  Systolic blood pressure

- tobacco:

  Cumulative tobacco consumption, in kg

- ldl:

  Low-density lipoprotein cholesterol

- adiposity:

  Adipose tissue concentration

- famhist:

  Family history of heart disease (1=Present, 0=Absent)

- typea:

  Score on test designed to measure type-A behavior

- obesity:

  Obesity

- alcohol:

  Current consumption of alcohol

- age:

  Age of subject

## Source

<https://web.stanford.edu/~hastie/ElemStatLearn/>

## References

- Hastie T, Tibshirani R, and Friedman J. (2001). *The Elements of
  Statistical Learning*. Springer.

- Rousseauw J, et al. (1983). Coronary risk factor screening in three
  rural communities. *South African Medical Journal*, **64**: 430-436.
