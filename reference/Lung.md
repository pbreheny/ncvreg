# VA lung cancer data set

Data from a randomised trial of two treatment regimens for lung cancer.
This is a standard survival analysis data set from the classic textbook
by Kalbfleisch and Prentice.

## Usage

``` r
Lung
```

## Format

A list of two objects: `y` and `X`

- y:

  A two column matrix (`Surv` object) containing the follow-up time (in
  days) and an indicator variable for whether the patient died while on
  the study or not.

- X:

  A matrix with 137 observations (rows) and 9 predictor variables
  (columns). The remainder of this list describes the columns of `X`

- trt:

  Treatment indicator (1=control group, 2=treatment group)

- karno:

  Karnofsky performance score (0=bad, 100=good)

- diagtime:

  Time from diagnosis to randomization (months)

- age:

  Age (years, at baseline)

- prior:

  Prior therapy (0=no, 1=yes)

- squamous:

  Indicator for whether the cancer type is squamous cell carcinoma
  (0=no, 1=yes)

- small:

  Indicator for whether the cancer type is small cell lung cancer (0=no,
  1=yes)

- adeno:

  Indicator for whether the cancer type is adenocarcinoma (0=no, 1=yes)

- large:

  Indicator for whether the cancer type is large cell carcinoma (0=no,
  1=yes)

## Source

<https://cran.r-project.org/package=survival>

## References

- Kalbfleisch D and Prentice RL (1980), *The Statistical Analysis of
  Failure Time Data*. Wiley, New York.

## See also

[`ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/ncvsurv.md)
