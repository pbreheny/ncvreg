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

## Passing in CV object --------------------------------------------------------
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- confidence_intervals(cv_fit)
run_tests(res)
check_pen(res, "lasso")
check_biased(res)

## Pass in ncvreg object -------------------------------------------------------
fit <- ncvreg(X, y, penalty = "SCAD")
expect_message({res <- confidence_intervals(fit)}, strict = TRUE)
run_tests(res)
check_pen(res, "SCAD")
check_biased(res)

## Pass in ncvreg object and X -------------------------------------------------
Xstd <- std(X)
fit <- ncvreg(Xstd, y, returnX = FALSE)
expect_message({res <- confidence_intervals(fit, X = Xstd)}, strict = TRUE)
run_tests(res)
check_biased(res)

## Expected errors -------------------------------------------------------------

## Pass in CV object and X and y, expect error as passed X is made different
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
expect_error(confidence_intervals(cv_fit, X = matrix(rnorm(500), 50, 10)))

## Check passing in non-standardized X (no error)
res <- confidence_intervals(cv_fit, X = X)
run_tests(res)

## And standardized X
res <- confidence_intervals(cv_fit, X = ncvreg::std(X))
run_tests(res)

## Pass in CV object with no X (expect error)
cv_fit <- cv.ncvreg(X, y, penalty = "lasso", returnX = FALSE)
expect_error(res <- confidence_intervals(cv_fit))

## Now supply X (don't expect error)
cv_fit <- cv.ncvreg(X, y, penalty = "lasso", returnX = FALSE)
res <- confidence_intervals(cv_fit, X = X)
run_tests(res)

## Pass some other object to fit
expect_error(confidence_intervals(lm(y ~ X)))

## Alternate Penalties ---------------------------------------------------------
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- confidence_intervals(cv_fit)
run_tests(res)
check_biased(res)
check_pen(res, "lasso")

cv_fit <- cv.ncvreg(X, y, penalty = "SCAD")
res <- confidence_intervals(cv_fit)
run_tests(res)
check_biased(res)
check_pen(res, "SCAD")

cv_fit <- cv.ncvreg(X, y, penalty = "SCAD", alpha = 0.5)
res <- confidence_intervals(cv_fit)
run_tests(res)
check_biased(res)
check_pen(res, "SCAD")
expect_equal(res$alpha[1], 0.5)

cv_fit <- cv.ncvreg(X, y, gamma = 4, alpha =0.7)
res <- confidence_intervals(cv_fit)
run_tests(res)
check_biased(res)
check_pen(res, "MCP")
expect_equal(res$gamma[1], 4)
expect_equal(res$alpha[1], 0.7)

## Alternate families
eta <- X %*% beta
y_bin <- rbinom(n = 50, size = 1, prob = exp(eta) / (1+exp(eta)))
cv_fit <- cv.ncvreg(X, y_bin, family = "binomial")
confidence_intervals(cv_fit)

y_pois <- rpois(50, exp(eta))
cv_fit <- cv.ncvreg(X, y_pois, family = "poisson", penalty = "MCP", alpha = 0.7)
confidence_intervals(cv_fit, adjust_projection = TRUE)

## LQA   -----------------------------------------------------------------------
## Passing in CV object
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- confidence_intervals(cv_fit, adjust_projection = TRUE)
run_tests(res)
check_biased(res)

## Pass in ncvreg object 
fit <- ncvreg(X, y, penalty = "SCAD")
res <- confidence_intervals(fit, adjust_projection = TRUE)
run_tests(res)
check_biased(res)

## Relaxed   -------------------------------------------------------------------
## Passing in CV object
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- confidence_intervals(cv_fit, relaxed = TRUE)
run_tests(res)

## Pass in ncvreg object 
fit <- ncvreg(X, y, penalty = "SCAD")
res <- confidence_intervals(fit, relaxed = TRUE)
run_tests(res)

## Pass in ncvreg object with poisson outcome
fit <- ncvreg(X, y_pois, penalty = "SCAD", family = "poisson")
res <- confidence_intervals(fit, relaxed = TRUE)
run_tests(res)

## Debiased -------------------------------------------------------------------
## Passing in CV object
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
res <- confidence_intervals(cv_fit, relaxed = TRUE, posterior = FALSE)
run_tests(res)
check_debiased(res)

## Pass in ncvreg object 
fit <- ncvreg(X, y, penalty = "SCAD")
res <- confidence_intervals(
  fit, relaxed = TRUE, adjust_projection = FALSE, posterior = FALSE
)
run_tests(res)
check_debiased(res)

## Lambda specification outside of range
fit <- ncvreg(X, y, penalty = "lasso")
lambda_seq <- fit$lambda
expect_error({
  confidence_intervals(fit, lambda = min(lambda_seq)*.5)
})
expect_error({
  confidence_intervals(fit, lambda = max(lambda_seq)*1.5)
})

## Run for null model
res <- confidence_intervals(fit, lambda = max(lambda_seq))

## Misc Checks
expect_error(confidence_intervals(ncvreg(X, y, penalty.factor = c(0, rep(1, 9)))))

## Examples
# Linear regression (SCAD-Net penalty, PIPE intervals, pass ncvreg object)
fit <- ncvreg(Prostate$X, Prostate$y, penalty = "SCAD", alpha = 0.9)
confidence_intervals(fit)

# Logistic regression (lasso penalty, LQA intervals, pass cv.ncvreg object)
data(Heart)
cv_fit <- cv.ncvreg(Heart$X, Heart$y, family="binomial", penalty = "lasso")
confidence_intervals(cv_fit, adjust_projection = TRUE) |> head()

## Singular issue warning
X[,2] <- 1
expect_error(ncvreg(X, y) |> confidence_intervals())
