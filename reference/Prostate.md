# Factors associated with prostate specific antigen

Data from a study by by Stamey et al. (1989) to examine the association
between prostate specific antigen (PSA) and several clinical measures
that are potentially associated with PSA in men who were about to
receive a radical prostatectomy.

## Usage

``` r
Prostate
```

## Format

A list of two objects: `y` and `X`

- y:

  Log PSA

- X:

  A matrix with 97 instances (rows) and 8 predictor variables (columns).
  The remainder of this list describes the columns of `X`

- lcavol:

  Log cancer volume

- lweight:

  Log prostate weight

- age:

  The man's age (years)

- lbph:

  Log of the amount of benign hyperplasia

- svi:

  Seminal vesicle invasion (1=Yes, 0=No)

- lcp:

  Log of capsular penetration

- gleason:

  Gleason score

- pgg45:

  Percent of Gleason scores 4 or 5

## Source

<https://web.stanford.edu/~hastie/ElemStatLearn/>

## References

- Hastie T, Tibshirani R, and Friedman J. (2001). *The Elements of
  Statistical Learning*. Springer.

- Stamey T, et al. (1989). Prostate specific antigen in the diagnosis
  and treatment of adenocarcinoma of the prostate. II. Radical
  prostatectomy treated patients. *Journal of Urology*, **16**:
  1076-1083.
