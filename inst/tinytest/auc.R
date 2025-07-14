if (interactive()) library(tinytest)
library(survival, quietly = TRUE)

# The quantities here are not exactly the same thing (one is in-sample, the
# other is out-of-sample), but they should be close
n <- 500
x <- matrix(rnorm(n * 10), n, 10)
y <- Surv(rexp(n, rate = exp(x[, 1])), rbinom(n, 1, prob = 0.8))
cvfit <- cv.ncvsurv(x, y, lambda.min = 0, returnY = TRUE)
a <- AUC(cvfit)
a[length(a)]
fit <- coxph(y ~ x)
b <- concordancefit(y, -predict(fit))$concordance
expect_equivalent(a[length(a)], b, tol = 0.02)
