rm(list = ls())
library(tinytest)

## Works
X <- matrix(rnorm(500), 50, 10)
beta <- c(1, rep(0, 9))
y <- X[,1] + rnorm(50)

run_tests <- function(res) {
  expect_inherits(res, "data.frame")
  expect_false(any(is.na(res[,c("lower", "upper")])))  
  expect_true(all(apply(res[,c("lower", "upper")], c(1, 2), is.finite)))  
}
check_pen <- function(res, pen = "MCP") {
  expect_equivalent(res$penalty[1], pen)
}
check_biased <- function(res) {
  expect_true(all(res$lower <= res$coef & res$coef <= res$upper))
}
check_debiased <- function(res) {
  expect_true(all(res$lower <= res$estimate & res$estimate <= res$upper))
}

## Passing data directly -------------------------------------------------------
expect_message({res <- pipe(X, y)}, strict = TRUE)
run_tests(res)
check_pen(res)
check_biased(res)

## Passing in CV object --------------------------------------------------------
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- pipe(fit = cv_fit)
run_tests(res)
check_pen(res, "lasso")
check_biased(res)

## Pass in ncvreg object -------------------------------------------------------
fit <- ncvreg(X, y, penalty = "SCAD")
expect_message({res <- pipe(fit = fit)}, strict = TRUE)
run_tests(res)
check_pen(res, "SCAD")
check_biased(res)

## Expected errors -------------------------------------------------------------

## Only pass X
expect_error(pipe(X))

## Pass in CV object and X and y, expect error as passed X is made different
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
expect_error(pipe(matrix(rnorm(500), 50, 10), y, fit = cv_fit))

## Check passing in non-standardized X (no error)
res <- pipe(X, y, fit = cv_fit)
run_tests(res)

## And standardized X
res <- pipe(ncvreg::std(X), y, fit = cv_fit)
run_tests(res)

## Check error with different y
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
expect_error(pipe(y = X[,1] + rnorm(50), fit = cv_fit))

## Pass in CV object with no X (expect error)
cv_fit <- cv.ncvreg(X, y, penalty = "lasso", returnX = FALSE)
expect_error(res <- pipe(fit = cv_fit))

## Now supply X (don't expect error)
cv_fit <- cv.ncvreg(X, y, penalty = "lasso", returnX = FALSE)
res <- pipe(X = X, fit = cv_fit)
run_tests(res)

## Alternate Penalties ---------------------------------------------------------
res <- pipe(X, y, penalty = "lasso")
run_tests(res)
check_biased(res)
check_pen(res, "lasso")

res <- pipe(X, y, penalty = "SCAD")
run_tests(res)
check_biased(res)
check_pen(res, "SCAD")

res <- pipe(X, y, penalty = "SCAD", alpha = 0.5)
run_tests(res)
check_biased(res)
check_pen(res, "SCAD")
expect_equal(res$alpha[1], 0.5)

res <- pipe(fit = cv.ncvreg(X, y, gamma = 4, alpha = 0.7))
run_tests(res)
check_biased(res)
check_pen(res, "MCP")
expect_equal(res$gamma[1], 4)
expect_equal(res$alpha[1], 0.7)

## Alternate families
eta <- X %*% beta
y_bin <- rbinom(n = 50, size = 1, prob = exp(eta) / (1+exp(eta)))
cv_fit <- cv.ncvreg(X, y_bin, family = "binomial")
pipe(fit = cv_fit)

y_pois <- rpois(50, exp(eta))
cv_fit <- cv.ncvreg(X, y_pois, family = "poisson", penalty = "lasso")
pipe(fit = cv_fit)

## LQA   -----------------------------------------------------------------------
## Passing data directly
expect_message({res <- pipe(X, y, adjust_projection = TRUE)}, strict = TRUE)
run_tests(res)
check_biased(res)

## Passing in CV object
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- pipe(fit = cv_fit, adjust_projection = TRUE)
run_tests(res)
check_biased(res)

## Pass in ncvreg object 
fit <- ncvreg(X, y, penalty = "SCAD")
expect_message({res <- pipe(fit = fit, adjust_projection = TRUE)}, strict = TRUE)
run_tests(res)
check_biased(res)

## Relaxed   -------------------------------------------------------------------
## Passing data directly
expect_message({res <- pipe(X, y, relaxed = TRUE)}, strict = TRUE)
run_tests(res)

## Passing in CV object
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- pipe(fit = cv_fit, relaxed = TRUE)
run_tests(res)

## Pass in ncvreg object 
fit <- ncvreg(X, y, penalty = "SCAD")
expect_message({res <- pipe(fit = fit, relaxed = TRUE)}, strict = TRUE)
run_tests(res)

## Debiased -------------------------------------------------------------------
## Passing data directly
expect_message({res <- pipe(X, y, posterior = FALSE, penalty = "lasso")}, strict = TRUE)
run_tests(res)
check_debiased(res)

## Passing in CV object
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- pipe(fit = cv_fit, relaxed = TRUE, posterior = FALSE)
run_tests(res)
check_debiased(res)

## Pass in ncvreg object 
fit <- ncvreg(X, y, penalty = "SCAD")
expect_message({res <- pipe(fit = fit, relaxed = TRUE, adjust_projection = FALSE, posterior = FALSE)}, strict = TRUE)
run_tests(res)
check_debiased(res)

## Lambda specification outside of range
lambda_seq <- ncvreg::ncvreg(X, y, penalty = "lasso")$lambda
expect_error({
  pipe(X, y, lambda = min(lambda_seq)*.5)
})
expect_error({
  pipe(X, y, lambda = max(lambda_seq)*1.5)
})

## Coercion 
expect_error(pipe(list(X), y))
res <- pipe(as.data.frame(X), y) ## A situation we would expect to work
run_tests(res)
res <- pipe(matrix(as.integer(X), 50, 10), y)
run_tests(res)

## Misc Checks
expect_error(pipe(fit = ncvreg(X, y, penalty.factor = c(0, rep(1, 9)))))

## Examples
data(Prostate)
X <- Prostate$X
y <- Prostate$y
pipe(X, y)

