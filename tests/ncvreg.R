require(ncvreg)
set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001, check.attributes=FALSE)}

#### Linear regression ####

# Works
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)
fit <- ncvreg(X, y, lambda.min=0, penalty="MCP")

# Equals MLE when lam=0
fit.mle <- lm(y~X)
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# logLik
stopifnot(equal(logLik(fit)[100], logLik(fit.mle)[1]))
stopifnot(equal(AIC(fit)[100], AIC(fit.mle)))

# Different penalties
fit <- ncvreg(X, y, lambda.min=0, penalty="SCAD")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))
fit <- ncvreg(X, y, lambda.min=0, penalty="lasso")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

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

#### Logistic regression ####

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

# Equals MLE when lam=0
fit.mle <- glm(y~X, family='binomial')
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# logLik
stopifnot(equal(logLik(fit)[100], logLik(fit.mle)[1]))
stopifnot(equal(AIC(fit)[100], AIC(fit.mle)))

# y logical
fit <- ncvreg(X, y==1, lambda.min=0, family='binomial')

#### Poisson regression ####

fit <- ncvreg(X, y, lambda.min=0, family='poisson')

# Equals MLE when lam=0
fit.mle <- glm(y~X, family='poisson')
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# logLik
stopifnot(equal(logLik(fit)[100], logLik(fit.mle)[1]))
stopifnot(equal(AIC(fit)[100], AIC(fit.mle)))

# Predict
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')
