suppressPackageStartupMessages(library(glmnet))

# ncvreg works for logistic regression
x <- matrix(rnorm(500), ncol = 10)
b <- rnorm(10)
y <- rnorm(x %*% b) > 0
beta <- glm(y ~ x, family = 'binomial')$coef
fit <- ncvreg(x, y, lambda = 1:0, family = 'binomial', penalty = 'SCAD', eps = 0.0001)
scad <- ncvreg(x, y, lambda = 1:0, family = 'binomial', penalty = 'SCAD', eps = 0.0001) |>
  coef(which = 2)
mcp <- ncvreg(x, y, lambda = 1:0, family = 'binomial', penalty = 'MCP', eps = 0.0001) |>
  coef(which = 2)
expect_equivalent(scad, beta, tolerance = 0.01)
expect_equivalent(mcp, beta, tolerance = 0.01)

# ncvreg reproduces lasso: binomial
nlasso <- coef(fit <- ncvreg(x, y, family = 'binomial', penalty = 'lasso'))
plot(fit, log = TRUE)
glasso <- glmnet(x, y, family = 'binomial', lambda = fit$lambda) |>
  coef() |>
  as.matrix()
plot(fit, 'lambda')
expect_equivalent(nlasso,  glasso, tolerance = 0.01)

# logLik() is correct: binomial
fit_mle <- glm(y ~ x, family = 'binomial')
fit <- ncvreg(x, y, lambda.min = 0, family = 'binomial')
expect_equivalent(logLik(fit)[100],  logLik(fit_mle)[1], tol = 0.001)
expect_equivalent(AIC(fit)[100], AIC(fit_mle), tol = 0.001)

# ncvreg dependencies work: binomial

# Predict
p <- predict(fit, x, 'link')
p <- predict(fit, x, 'response')
p <- predict(fit, x, 'class')
p <- predict(fit, x, 'coef')
p <- predict(fit, x, 'vars')
p <- predict(fit, x, 'nvars')

# y logical
fit <- ncvreg(x, y == 1, lambda.min = 0, family = 'binomial')

# Summary
summary(fit, which = 10)
summary(fit, lam = 0.05)

# cv.ncvreg() options work for binomial
n <- 200
p <- 10
x <- matrix(rnorm(n * p), ncol = p)
b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean = x %*% b, sd = 1) > 0.5

par(mfrow = c(2, 2))
cvfit <- cv.ncvreg(x, y, family = 'binomial')
plot(cvfit, type = 'all')
summary(cvfit)
predict(cvfit, type = 'coefficients')
predict(cvfit, type = 'vars')
predict(cvfit, type = 'nvars')
head(predict(cvfit, x = x, 'link'))
head(predict(cvfit, x = x, 'response'))
head(predict(cvfit, x = x, 'class'))

b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean = x %*% b, sd = 5) > 0.5
cvfit <- cv.ncvreg(x, y, family = 'binomial')
plot(cvfit, type = 'all')

b <- rep(0, 10)
y <- rnorm(n, mean = x %*% b, sd = 5) > 0.5
cvfit <- cv.ncvreg(x, y, family = 'binomial')
plot(cvfit, type = 'all')
