# Package index

## Model fitting

- [`ncvfit()`](https://pbreheny.github.io/ncvreg/reference/ncvfit.md) :
  Direct interface for nonconvex penalized regression (non-pathwise)
- [`ncvreg()`](https://pbreheny.github.io/ncvreg/reference/ncvreg.md) :
  Fit an MCP- or SCAD-penalized regression path
- [`ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/ncvsurv.md)
  : Fit an MCP- or SCAD-penalized survival model
- [`std()`](https://pbreheny.github.io/ncvreg/reference/std.md) :
  Standardizes a design matrix

## Cross-validation

- [`assign_fold()`](https://pbreheny.github.io/ncvreg/reference/assign_fold.md)
  : Assign folds for cross-validation
- [`cv.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.md)
  [`cv.ncvsurv()`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.md)
  : Cross-validation for ncvreg/ncvsurv
- [`plot(`*`<cv.ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/plot.cv.ncvreg.md)
  : Plots the cross-validation curve from a cv.ncvreg object
- [`summary(`*`<cv.ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/summary.cv.ncvreg.md)
  [`print(`*`<summary.cv.ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/summary.cv.ncvreg.md)
  : Summarizing cross-validation-based inference
- [`AUC(`*`<cv.ncvsurv>`*`)`](https://pbreheny.github.io/ncvreg/reference/AUC.md)
  : AUC for cv.ncvsurv objects

## Plotting and extracting model features

- [`logLik(`*`<ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/logLik.ncvreg.md)
  [`logLik(`*`<ncvsurv>`*`)`](https://pbreheny.github.io/ncvreg/reference/logLik.ncvreg.md)
  : Extract Log-Likelihood

- [`predict(`*`<cv.ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/predict.ncvreg.md)
  [`coef(`*`<cv.ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/predict.ncvreg.md)
  [`predict(`*`<cv.ncvsurv>`*`)`](https://pbreheny.github.io/ncvreg/reference/predict.ncvreg.md)
  [`predict(`*`<ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/predict.ncvreg.md)
  [`coef(`*`<ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/predict.ncvreg.md)
  : Model predictions based on a fitted ncvreg object.

- [`predict(`*`<ncvsurv>`*`)`](https://pbreheny.github.io/ncvreg/reference/predict.ncvsurv.md)
  :

  Model predictions based on a fitted `ncvsurv` object.

- [`plot(`*`<ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/plot.ncvreg.md)
  : Plot coefficients from a ncvreg object

- [`plot(`*`<ncvsurv.func>`*`)`](https://pbreheny.github.io/ncvreg/reference/plot.ncvsurv.func.md)
  : Plot survival curve for ncvsurv model

- [`residuals(`*`<ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/residuals.ncvreg.md)
  : Extract residuals from a ncvreg or ncvsurv fit

## Inference

- [`local_mfdr()`](https://pbreheny.github.io/ncvreg/reference/local_mfdr.md)
  : Estimate local mFDR for all features
- [`mfdr()`](https://pbreheny.github.io/ncvreg/reference/mfdr.md) :
  Marginal false discovery rates
- [`plot(`*`<mfdr>`*`)`](https://pbreheny.github.io/ncvreg/reference/plot.mfdr.md)
  : Plot marginal false discovery rate curves
- [`perm.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/perm.ncvreg.md)
  : Permutation fitting for ncvreg
- [`permres()`](https://pbreheny.github.io/ncvreg/reference/permres.md)
  : Permute residuals for a fitted ncvreg model
- [`summary(`*`<ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/summary.ncvreg.md)
  [`print(`*`<summary.ncvreg>`*`)`](https://pbreheny.github.io/ncvreg/reference/summary.ncvreg.md)
  : Summary method for ncvreg objects
- [`intervals()`](https://pbreheny.github.io/ncvreg/reference/intervals.md)
  : Projection base test statistics and intervals

## Data sets

- [`Heart`](https://pbreheny.github.io/ncvreg/reference/Heart.md) : Risk
  factors associated with heart disease
- [`Lung`](https://pbreheny.github.io/ncvreg/reference/Lung.md) : VA
  lung cancer data set
- [`Prostate`](https://pbreheny.github.io/ncvreg/reference/Prostate.md)
  : Factors associated with prostate specific antigen
