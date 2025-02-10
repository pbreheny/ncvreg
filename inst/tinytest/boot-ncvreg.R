rm(list = ls())
library(tinytest)

## Works
X <- matrix(rnorm(500), 50, 10)
y <- X[,1] + rnorm(50)

run_tests <- function(boot_res) {
  expect_inherits(boot_res, "list")
  expect_inherits(boot_res$confidence_intervals, "data.frame")
  expect_equivalent(dim(boot_res$confidence_intervals), c(10, 3))
  expect_false(any(is.na(boot_res$confidence_intervals)))  
  expect_true(all(apply(boot_res$confidence_intervals[,c("lower", "upper")], c(1, 2), is.finite)))  
}

## Regular usage, also test returning boot draws
boot_res <- boot_ncvreg(X, y, return_boot = TRUE)
run_tests(boot_res)
expect_equal(dim(boot_res$boot_draws), c(1000, 10))

## Pass in CV object
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
boot_res <- boot_ncvreg(fit = cv_fit, returnCV = TRUE)
run_tests(boot_res)
expect_equivalent(cv_fit, boot_res$cv.ncvreg)

## Pass in ncvreg object
fit <- ncvreg(X, y, penalty = "lasso")
boot_res <- boot_ncvreg(fit = fit)
run_tests(boot_res)

## Pass in X but not y or fit object, expect error
expect_error(boot_ncvreg(X))

## Check if seed seeting is working
seed_before <- .GlobalEnv$.Random.seed
boot_res <- boot_ncvreg(X, y, seed = 1234)
seed_after <- .GlobalEnv$.Random.seed
expect_identical(seed_before, seed_after)

## Pass in CV object and X and y, expect error as passed X is made different
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
expect_error(boot_ncvreg(matrix(rnorm(500), 50, 10), y, fit = cv_fit))

## Check passing in non-standardized X
boot_res <- boot_ncvreg(X, y, fit = cv_fit)
run_tests(boot_res)

## And standardized X
boot_res <- boot_ncvreg(ncvreg::std(X), y, fit = cv_fit)
run_tests(boot_res)

## Check error with different y
cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
expect_error(boot_ncvreg(y = X[,1] + rnorm(50), fit = cv_fit))

## Pass in CV object with no X (expect error)
cv_fit <- cv.ncvreg(X, y, penalty = "lasso", returnX = FALSE)
expect_error(boot_res <- boot_ncvreg(fit = cv_fit))

## Now supply X (don't expect error)
cv_fit <- cv.ncvreg(X, y, penalty = "lasso", returnX = FALSE)
boot_res <- boot_ncvreg(X = X, fit = cv_fit)
run_tests(boot_res)

## Alternate Penalties
boot_res <- boot_ncvreg(X, y, penalty = "MCP")
run_tests(boot_res)

boot_res <- boot_ncvreg(X, y, penalty = "SCAD")
run_tests(boot_res)

boot_res <- boot_ncvreg(X, y, penalty = "lasso", alpha = 0.5)
run_tests(boot_res)

## Lambda specification outside of range
lambda_seq <- ncvreg::ncvreg(X, y, penalty = "lasso")$lambda
expect_message({
  boot_ncvreg(X, y, lambda = min(lambda_seq))
}, strict = TRUE)
expect_message({
  boot_ncvreg(X, y, lambda = max(lambda_seq))
}, strict = TRUE)

## Coercion 
expect_error(boot_ncvreg(list(X), y))
boot_res <- boot_ncvreg(as.data.frame(X), y) ## A situation we would expect to work
run_tests(boot_res)
boot_res <- boot_ncvreg(matrix(as.integer(X), 50, 10), y)
run_tests(boot_res)

## Misc Checks
expect_warning(boot_ncvreg(X, fit = cv_fit, nlambda = 10), strict = TRUE)
expect_error(boot_ncvreg(X, y, penalty.factor = c(0, rep(1, 9))))
expect_error(boot_ncvreg(X, y, family = "poisson"))
expect_error({
  y_bin <- exp(y) / (1+exp(y)) > 0.5
  cv_fit <- cv.ncvreg(X, y_bin, family = "binomial")
  boot_ncvreg(fit = cv_fit)
})
expect_message(boot_ncvreg(X, y, verbose = TRUE), strict = TRUE)
expect_message(boot_ncvreg(X, y, verbose = TRUE, sigma2 = 1, lambda = 0.1), strict = TRUE)
expect_message(boot_ncvreg(X, y, returnX = TRUE, convex = TRUE), strict = TRUE)

## Parallelization support
cl <- parallel::makeCluster(4)
runtime_orig <- system.time({
  boot_res <- boot_ncvreg(X, y, seed = 123445)
})
runtime_par <- system.time({
  boot_res_par <- boot_ncvreg(X, y, cluster = cl)
})
parallel::stopCluster(cl)
run_tests(boot_res_par)
expect_true(runtime_par[[2]] < runtime_orig[[2]])

cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
expect_error(boot_ncvreg(fit = cv_fit, cluster = NA)) ## Passing something other than a cluster

## Nonsingular issue
X[,2] <- 1
expect_warning({
  boot_res <- boot_ncvreg(X, y)
}, strict = TRUE)
expect_true(all(is.na(boot_res$confidence_intervals[2,c("lower", "upper")])))

## Nonsingular issue with cv.ncvreg object passed
expect_warning({
  cv_fit <- cv.ncvreg(X, y, penalty = "lasso")
  boot_res <- boot_ncvreg(fit = cv_fit)
}, strict = TRUE)
expect_true(all(is.na(boot_res$confidence_intervals[2,c("lower", "upper")])))

## Examples
data(Prostate)
X <- Prostate$X
y <- Prostate$y
boot_ncvreg(X, y, max.iter = 100)
boot_ncvreg(fit = cv.ncvreg(X, y, penalty = "lasso"))

ncvreg(X, y, penalty = "lasso", max.iter = 100)
