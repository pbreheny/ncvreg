if (interactive()) library(tinytest)
suppressPackageStartupMessages(library(parallel))

## lasso
X <- matrix(rnorm(500), 50, 10)
y <- X[,1] + rnorm(50)
boot <- boot_ncvreg(X, y)

## Change bootstrap sample numbers
boot <- boot_ncvreg(X, y, nboot = 100)

## Set seed
boot_seed <- boot_ncvreg(X, y, seed = 1234, nboot = 100)

## MCP
boot <- boot_ncvreg(X, y, penalty = "MCP", nboot = 100)

## SCAD
boot <- boot_ncvreg(X, y, penalty = "SCAD", nboot = 100)

## Elastic Net
boot <- boot_ncvreg(X, y, alpha = 0.5, nboot = 100)

## Supply cv.ncvreg object directly
set.seed(1234)
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
boot <- boot_ncvreg(cv_fit = cv_fit, nboot = 100)
tinytest::expect_equal(boot, boot_seed)

## Parallel
cl <- makeCluster(4)
boot <- boot_ncvreg(X, y, cluster = cl, nboot = 100)

## User supplied lambda
boot <- boot_ncvreg(X, y, nboot = 100, lambda = cv_fit$lambda.min)

## User supplied sigma2
boot <- boot_ncvreg(X, y, nboot = 100, sigma2 = 1)

## User supplied lambda and sigma2
boot <- boot_ncvreg(X, y, nboot = 100, lambda = cv_fit$lambda.min, sigma2 = 1)

## Return CV fit
boot <- boot_ncvreg(X, y, nboot = 100, returnCV = TRUE)
tinytest::expect_equivalent(class(boot$cv.ncvreg), "cv.ncvreg")

## Expected errors
tinytest::expect_error(boot_ncvreg(penalty = "lasso")) ## missing X and y
cv_fit_noX <- cv.ncvreg(X, y, returnX = FALSE, penalty = "lasso")
tinytest::expect_error(boot_ncvreg(cv_fit = cv_fit_noX)) ## No X in cv_fit
y_binom <- rbinom(n = 50, 1, prob = 0.5)
cv_fit_binom <- cv.ncvreg(X, y_binom, family = "binomial", penalty = "lasso")
tinytest::expect_error(boot_ncvreg(cv_fit = cv_fit_binom)) ## No X in cv_fit
tinytest::expect_error(boot_ncvreg(X = cv_fit, y)) ## X not coercable type
tinytest::expect_error(boot_ncvreg(X, y = cv_fit)) ## y not coercable type
tinytest::expect_error(boot_ncvreg(X, y, penalty.factor = 0.5)) ## y not coercable type
tinytest::expect_error(boot_ncvreg(cv_fit = cv_fit, lambda = 0.5*min(cv_fit$lambda))) ## Trying to estimate sigma2 for lambda outside of fit
tinytest::expect_error(boot_ncvreg(X, y, cluster = "cl")) ## Trying to estimate sigma2 for lambda outside of fit

## Expected warnings
tinytest::expect_warning(boot_ncvreg(X, y, family = "binomial")) ## Attempt to apply to non-Gaussian outcome
tinytest::expect_warning(boot_ncvreg(X, y, cv_fit)) ## Attempt to apply to non-Gaussian outcome
tinytest::expect_warning(boot_ncvreg(cv_fit = cv_fit, warn = TRUE)) ## Attempt to supply additional arguments

## Expect messages
tinytest::expect_message(boot_ncvreg(X, y, verbose = FALSE, nboot = 100, lambda = min(cv_fit$lambda) * 0.5)) ## lambda too smalle
tinytest::expect_message(boot_ncvreg(X, y, verbose = FALSE, nboot = 100, lambda = max(cv_fit$lambda) * 2)) ## lambda too large
tinytest::expect_message(boot_ncvreg(cv_fit = cv_fit, lambda = cv_fit$lambda.min, sigma2 = 1)) ## Ignore user supplied value
tinytest::expect_message(boot_ncvreg(cv_fit = cv_fit, lambda = cv_fit$lambda.min)) ## Estimating sigma2 at supplied value of lambda

## Coercion
X_int <- matrix(rnorm(500), 50, 10)
storage.mode(X_int) <- "integer"
boot <- boot_ncvreg(X_int, y, nboot = 100)

## Applied to prostate
data(Prostate)
tinytest::expect_message(boot_ncvreg(Prostate$X, Prostate$y, penalty = "lasso", verbose = FALSE, nboot = 100, lambda = min(cv_fit$lambda) * 0.5))
