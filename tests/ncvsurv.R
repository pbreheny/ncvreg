require(ncvreg)
require(survival)
set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001, check.attributes=FALSE)}

# Works
y <- Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*10), 50, 10)
fit <- ncvsurv(X, y, lambda.min=0)

# Equals MLE when lam=0
fit.mle <- coxph(y~X)
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# logLik
stopifnot(equal(logLik(fit)[100], logLik(fit.mle)[1]))
stopifnot(equal(AIC(fit)[100], AIC(fit.mle)))

# Other penalties
fit <- ncvsurv(X, y, lambda.min=0, penalty="SCAD")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))
fit <- ncvsurv(X, y, lambda.min=0, penalty="lasso")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# Predict
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')
S <- predict(fit, X, 'survival', lambda=0.1)
plot(S)
S <- predict(fit, X[1,], 'survival', lambda=0.1)
plot(S)
p <- predict(fit, X, 'median')

# Penalty factor
fit <- ncvsurv(X, y, penalty.factor=c(0:9))

# Coersion
fit <- ncvsurv(as.data.frame(X), y)

# Loss
eta <- predict(fit, X, 'link', lambda=0.1)
ncvreg:::loss.ncvsurv(y, eta)
