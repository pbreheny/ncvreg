# ncvreg 3.15.0
  * New: boot_ncvreg() function to obtain confidence intervals
  * New: assign_fold() function to assign folds for CV
  * Change: seed is no longer an argument to CV functions; use assign_fold()
    instead

# ncvreg 3.14.3
  * Internal: Now using R_Calloc for _R_USE_STRICT_R_HEADERS_ compatibility

# ncvreg 3.14.2
  * Documentation: Lots of formatting fixes to the documentation

# ncvreg 3.14.1
  * Fixed: cv.ncvreg(), cv.ncvsurv() no longer affect seed in global
    environment if seed is specified

# ncvreg 3.14.0
  * New: residuals() method
  * New: std() can now be applied to new data
  * New: summary.ncvreg() now offers sort option; fixes #13
  * Change: fir() deprecated
  * Change: local_mfdr() allows user to specify sigma; also uses CV if
    called with cv object
  * Fixed: Manual color palettes now recycled correctly; fixes #40;
    thank you to Logan Harris for pointing this out
  * Fixed: mfdr now works for Poisson
  * Documentation: Adding vignettes on other CV criteria, adaptive
    rescaling
  * Documentation: References reformatted, URLs updated, DOIs added
  * Internal: C code for binomial, poisson now unified under glm
    structure
  * Internal: Now using roxygen2 for all documentation

# ncvreg 3.13.0
  * New: Options 'xtx' and 'r' for ncvfit()
  * Internal: cv.ncvreg() now uses less memory (returnX off)
  * Internal: Better error handling if a matrix is supplied for y
  * Fixed: AUC() now compatible with survival 3.2.10

# ncvreg 3.12.0
  * New: ncvfit(), a raw API to the ncvreg solver with full control over
    standardization, etc.
  * Changed: ncvreg and ncvsurv now issue warning for non-pathwise usage
  * Internal: Now using tinytest for unit testing
  * Fixed: Memory leak in cox-dh; resolves #20

# ncvreg 3.11.2
  * New: std() now works on integer matrices and numeric vectors
  * Internal: Lots of internal changes for cleaner, more reliable code
  * New version numbering system

# ncvreg 3.11-1
  * Fixed: Leave-one-out cross-validation now works correctly for logistic
    regression
  * Documentation: Added documentation (online) for local mfdr
  * Documentation: Fixed some broken links and typos

# ncvreg 3.11-0
  * Change: returnX now turned on by default if X < 100 Mb (used to be 10 Mb)
  * Change: summary.ncvreg now based solely on local mfdr
  * Change: Loss functions now consistently defined as deviance for all types
    of models
  * Change: R^2 now consistently uses the Cox-Snell definition for all types
    of models
  * Change: cv.ncvreg and cv.ncvsurv now return fold assignments
  * Fixed: Can now pass fold assignments to cv.ncvsurv
  * Documentation: Lots of updates
  * Documentation: vignette now html (used to be pdf)
  * Documentation: pkgdown website

# ncvreg 3.10-0
  * New: summary.ncvreg and summary.ncvsurv now report tables of inference for
    each feature based on local mFDRs
  * New: Option to specify fold assignments in cv.ncvsurv
  * New: CVSE now calculated for Cox models, with option of quick or bootstrap
  * Change: returnX now turned on by default if X < 10 Mb
  * Change: cv.ncvsurv now balances censoring across fold assignments
  * Change: All data sets now follow Data$X, Data$y convention
  * Deprecated: cv.ind argument to cv.ncvreg is now called fold
  * Portability: Fixed C99 flag
  * Internal: Fixed & v && C issue

# ncvreg 3.9-1
  * Change: Poission now returns linear predictors, like other families
  * Internal: Changing PROTECT/UNPROTECT to conform to new coding standards

# ncvreg 3.9-0
  * Deprecated: fir() is now called mfdr()
  * Change: mfdr for Cox and logistic models no longer use the simplistic
    approximation of 3.7-0.  These calculations are much more accurate, but
    more computationally intensive, so these are carried out in C now.
  * Change: mfdr for Cox and logistic models requires the model matrix X now.
  * Internal: Registration of native routines
  * Fixed: std() wasn't matching up column names if one column got dropped

# ncvreg 3.8-0
  * Change: max.iter now based on total number of iterations for entire path
  * Fixed: Bug when fitting Cox model for single lambda
  * Fixed: std no longer drops dimnames

# ncvreg 3.7-1
  * Fixed: Various fixes for fir function
  * Fixed: Bug with high dimensional (p > n) Cox models

# ncvreg 3.7-0
  * New: fir extended to Cox and logistic regression
  * New: summary function for ncvreg and ncvsurv objects
  * Change: Convergence criterion now based on RMSD of linear predictors
  * Change: Additional options and improvements to plot.fir
  * Change: Better display of fir objects
  * Internal: Improved efficiency for Cox models (linear predictor calculation
    now occurs in C, not R)
  * Internal: Reorganized testing suite
  * Fixed: lamNames with single lambda passed
  * Fixed: loss wasn't being returned for gaussian if failure to converge
  * Fixed: perm.ncvreg would return NAs when models were saturated

# ncvreg 3.6-0
  * New: Exports std() function for standardizing a design matrix
  * Fixed: In predict.cv.ncvsurv
  * Documentation: Added 'quick start' vignette
  * Internal: Improved efficiency for cox models (avoids recalculating linear
    predictors)
  * Internal: Reorganized testing suite
  * Internal: 'survival' package now used for setupLambda in Cox models

# ncvreg 3.5-2
  * New: Added user interrupt checking
  * Fixed: In ncvsurv with integer penalty factors
  * Fixed: Rare numerical accuracy bug in cv fold assignments
  * Fixed: LOOCV bug introduced by bias-correction feature

# ncvreg 3.5-1
  * New: Compute bias correction for CV error; this is an experimental
    feature at this point and may change in the future
  * Internal: Replaced AUC function with more efficient version using
    survival package
  * Fixed: Penalty.factor for cv.ncvsurv when some columns may be degenerate

# ncvreg 3.5-0
  * New: Added function AUC() to calculate cross-validated AUC values
    for ncvsurv models.
  * New: Option to return fitted values from cross-validation folds
    (returnY=TRUE) for cv.ncvreg and cv.ncvsurv.
  * Change: New method for calculation of cross-validation loss in cv.ncvsurv.
  * Change: More accurate calculation for convexMin in the presence of
    unpenalized variables
  * Fixed: Factor-valued y with CV logistic regression
  * Internal: Substantial efficiency improvements throughout for Cox models.
    Coordinate descent redesigned to work in O(n) instead of O(n^2) operations,
    and R code redesigned at various points to avoid the creation of any n x n
    matrices when fitting and cross-validating Cox regression models.
  * Internal: Better double/int type checking for penalty.factor
  * Internal: Modifications to NAMESPACE for compatibility with R 3.3.

# ncvreg 3.4-0
  * New: Expanded predict function for Cox models.  predict.ncvsurv now
    estimates subject-specific survival functions and medians.
  * New: Plot method for survival curves.
  * New: Option in perm.ncvreg to permute residuals for linear
    regression
  * New: permres function to estimate false inclusion rates based on
    residuals at a specific value of lambda
  * New: Some support for factors in X, y.  It is still recommended that
    users convert X to a numeric matrix prior to fitting in order to ensure that
    predict() methods work properly, but ncvreg will now allow you to pass a
    data frame with factors and handle things appropriately.
  * Fixed: In predict.ncvsurv, when applied to models with saturation issues.
  * Fixed: Small memory leak in ncvsurv.

# ncvreg 3.3-0
  * New: Support for fitting survival models added (ncvsurv), along
    with predict, plot, and cv.ncvsurv support functions.  Currently, Cox models
    are the only type of survival model implemented.
  * New: Parallelization support for cv.ncvreg (with help from Grant
    Brown)
  * Fixed: In cv.ncvreg, when attempting to use leave-one-out cross-validation
    (thank you to Cajo ter Braak for pointing this out)
  * Removed: ncvreg_fit; it may return in a future version of the package.

# ncvreg 3.2-0
  * New: Automatically coerces X to matrix and y to numeric if possible
  * New: Made ncvreg_fit more user-friendly: user no longer has to specify
    lambda, works with coef, predict, plot, etc.
  * Changed: Modified order of arguments for predict so that 'type' comes
    before 'lambda' and 'which'
  * Fixed: Bug in convexMin when used with penalty.factor option
  * Internal: Updated algorithm to 'hybrid' strong/active cycling

# ncvreg 3.1-0
  * New: Added support for Poisson regression
  * Fixed: Bug in ncvreg_fit that could arise when fitting a model without an
    intercept
  * Fixed: Bug in cv.ncvreg with univariate regression (thank you to Diego
    Franco Saldana for pointing this out)

# ncvreg 3.0-0
  * New: Added fir, perm.ncvreg, and plot.fir functions for the purposes of
    estimating and displaying false inclusion rates; these are likely to evolve
    over the next few months
  * Fixed: Bug in cv.ncvreg for user-specified lambda sequence
  * Internal: Revised algorithms to incorporate targeted cycling based on strong
    rules
  * Internal: Moved standardization to C
  * Internal: Moved calculation of lambda sequence to C
  * Internal: As a result of the above three changes, ncvreg now runs much
    faster for large p

# ncvreg 2.7-0
  * New: "vars" and "nvars" options to predict function.
  * Changed: Modified look of summary(cvfit) output.
  * Internal: Modified details of .Call interface.

# ncvreg 2.6-0
  * New: Introduction of function ncvreg_fit for programmers who want to access
    the internal C routines of ncvreg, bypassing internal standardization and
    processing
  * New: Added vertical.line and col options to plot.cv.ncvreg
  * Fixed: Bug in axis annotations with plot.cv.ncvreg when model is saturated
  * Fixed: Deviance calculation; would return NaN if fitted probabilities of 0
    or 1 occurred for binomial outcomes
  * Fixed: NAMESPACE for coef.cv.ncvreg and predict.cv.ncvreg
  * Internal: .Call now used instead of .C

# ncvreg 2.5-0
  * New: Options in plot.cv.ncvreg to plot estimates of r-squared,
    signal-to-noise ratio, scale parameter, and prediction error in addition to
    cross-validation error (deviance)
  * New: Summary method for cv.ncvreg which displays the above information at
    lambda.min, the value of lambda minimizing the cross-validation error
  * Fixed: Bug in cv.ncvreg with user-defined lambda values.

# ncvreg 2.4-0
  * New: penalty.factor option
  * New: coef and predict methods now accept lambda as argument
  * New: logLik method (which in turn allows AIC/BIC)
  * Changed: cv.grpreg now returns full data fit as well as CV errors
  * Fixed: Error in definition/calculation of cross-validation error and
    standard error
  * Fixed: Bug that arose if lambda was scalar (instead of a vector)
  * Fixed: Bug in cv.ncvreg for linear regression -- cross-validation was being
    carried out deterministically (Thank you to Brenton Kenkel for pointing this
    out)
  * Fixed: Intercept for logistic regression was not being calculated for
    lamda=0
  * Internal: standardization more efficient
  * Internal: cdfit_ now returns loss (RSS for gaussian, deviance for binomial)

# ncvreg 2.3-2
  * Documentation: Fixed formatting error in citation.
  
# ncvreg 2.3-1
  * Changed: plot.ncvreg: Made the passing of arguments for plot.ncvreg more
    flexible, so that user can pass options concerning both the plot and the
    lines
  * Changed: plot.ncvreg: Changed some of the default settings with respect to
    color (hcl instead of hsv) and line width

# ncvreg 2.3
  * Documentation: Updated documentation for cv.ncvreg.Rd, which no longer
    agreed with the function usage (this was an oversight in the release of
    version 2.2)

# ncvreg 2.2
  * New: plot.cv.ncvreg for plotting cv.ncvreg objects
  * Changed: Divorced cross-validation from fitting in cv.ncvreg.  From a user
    perspective, this increases flexibility, although obtaining the model with
    CV-chosen regularization parameter now requires two calls (to ncvreg and
    cv.ncvreg).  The functions, however, are logically separate and involve
    entirely separate methods.
