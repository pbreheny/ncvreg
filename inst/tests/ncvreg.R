set.seed(1)

# Linear regression -------------------------------------------------------

# Works
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)
fit <- ncvreg(X, y, lambda.min=0, penalty="MCP")

.test="Equals MLE when lam=0"
fit.mle <- lm(y~X)
check(coef(fit)[,100], coef(fit.mle), tol=0.001)
fit <- ncvreg(X, y, lambda.min=0, penalty="SCAD")
check(coef(fit)[,100], coef(fit.mle), tol=0.001)
fit <- ncvreg(X, y, lambda.min=0, penalty="lasso")
check(coef(fit)[,100], coef(fit.mle), tol=0.001)

.test="logLik is correct"
check(logLik(fit)[100], logLik(fit.mle)[1])
check(AIC(fit)[100], AIC(fit.mle), tol=0.001)

# Predict
p <- predict(fit, X, 'link', lambda=0.1)
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')

# Integers
X <- matrix(rpois(500, 1), 50, 10)
y <- rpois(50, 1)
fit <- ncvreg(X, y)

# Data frame
fit <- ncvreg(as.data.frame(X), y)

# Penalty factor
fit <- ncvreg(as.data.frame(X), y, penalty.factor=c(0:9))

# User lambdas
fit <- ncvreg(as.data.frame(X), y, lambda=c(1, 0.1, 0.01))

# ReturnX
fit <- ncvreg(as.data.frame(X), y, returnX=TRUE)

# Constant columns
fit <- ncvreg(cbind(5, X), y)

# Plot
plot(fit)
plot(fit, log.l=TRUE)

# Summary
summary(fit, which=10)
summary(fit, lam=0.05)

# Logistic regression -----------------------------------------------------

# Works
y <- rbinom(50, 1, 0.5)
fit <- ncvreg(X, y, lambda.min=0, family='binomial')

# Predict
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'class')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')

.test="Equals MLE (logistic)"
fit.mle <- glm(y~X, family='binomial')
check(coef(fit)[,100], coef(fit.mle), tol=0.001)

.test="logLik is correct (logistic)"
check(logLik(fit)[100], logLik(fit.mle)[1], tol=0.001)
check(AIC(fit)[100], AIC(fit.mle), tol=0.001)

# y logical
fit <- ncvreg(X, y==1, lambda.min=0, family='binomial')

# Summary
summary(fit, which=10)
summary(fit, lam=0.05)

# Poisson regression ------------------------------------------------------

fit <- ncvreg(X, y, lambda.min=0, family='poisson')

.test="Equals MLE (Poisson)"
fit.mle <- glm(y~X, family='poisson')
check(coef(fit)[,100], coef(fit.mle), tol=0.001)

.test="logLik is correct (Poisson)"
check(logLik(fit)[100], logLik(fit.mle)[1], tol=0.001)
check(AIC(fit)[100], AIC(fit.mle), tol=0.001)

# Predict
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')
