# Single lambda warning
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)
expect_warning(ncvreg(X, y, lambda = 0.1))
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
expect_warning(ncvsurv(X, y, lambda = 0.1))

# mfdr
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n) > 0
fit <- lm(y ~ X)
expect_error(mfdr(fit))
fit <- ncvreg(X, y, family = "binomial", returnX = FALSE)
expect_error(mfdr(fit))

# ncvreg
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)
expect_error(ncvreg(rnorm, y), "must be a matrix")
expect_error(ncvreg(matrix(LETTERS[1:10], 10, 1), 1:10), "numeric matrix")
expect_error(ncvreg(X, "a"), "must be numeric")
expect_error(ncvreg(X, y, penalty = "MCP", gamma = 0.5), "greater than 1")
expect_error(ncvreg(X, y, penalty = "SCAD", gamma = 1.5), "greater than 2")
expect_error(ncvreg(X, y, nlambda = 1), "at least 2")
expect_error(ncvreg(X, y, alpha = 0), "greater than 0")
expect_error(ncvreg(X, c(NA, y[-1])), "missing data")
expect_error(ncvreg(X, y, penalty.factor = 1:floor(p / 2)), "penalty.factor and X do not match")
expect_error(ncvreg(X, y, family = "binomial"), "non-binary")

# ncvfit
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)
expect_error(ncvfit(rnorm, y), "must be a matrix")
expect_error(ncvfit(matrix(LETTERS[1:10], 10, 1), 1:10, lambda = 0.1), "numeric matrix")
expect_error(ncvfit(X, "a"), "must be numeric")
expect_error(ncvfit(X, y, penalty = "MCP", gamma = 0.5), "greater than 1")
expect_error(ncvfit(X, y, penalty = "SCAD", gamma = 1.5), "greater than 2")
expect_error(ncvfit(X, y, alpha = 0), "greater than 0")
expect_error(ncvfit(X, c(NA, y[-1])), "missing data")
expect_error(ncvfit(X, y, penalty.factor = 1:floor(p / 2)), "penalty.factor and X do not match")

# permres
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)
fit <- ncvreg(X, y)
expect_error(permres(fit, lambda = "a"), "must be numeric")
fit_nox <- ncvreg(X, y, returnX = FALSE)
expect_error(permres(fit_nox, lambda = 0.1), "returnX")
