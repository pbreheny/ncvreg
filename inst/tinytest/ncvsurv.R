library(survival, quietly=TRUE)
if (interactive()) library(tinytest)

# ncvsurv works for simple cox regression, no censoring -------------------

n <- 10
X <- matrix(runif(n), n, 1)
y <- Surv(sort(rexp(n, X), decreasing=TRUE), rep(1,n))
X <- X - mean(X)
beta <- coef(coxph(y ~ X))
eta <- X%*%beta
w <- exp(eta)/cumsum(exp(eta))
sum(X * (1-w))
scad <- coef(ncvsurv(X, y, lambda.min=0, penalty="SCAD", eps=.0001), lambda=0)
mcp <- coef(ncvsurv(X, y, lambda.min=0, penalty="MCP", eps=.0001), lambda=0)
expect_equivalent(scad, beta, tolerance=.01)
expect_equivalent(mcp, beta, tolerance=.01)


# ncvsurv works for simple cox regression ---------------------------------

n <- 30
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
y[which.max(y[,1]), 2] <- 0
y[which.min(y[,1]), 2] <- 0
X <- matrix(rnorm(n), n, 1)
beta <- coef(coxph(y ~ X))
scad <- coef(ncvsurv(X, y, lambda.min=0, penalty="SCAD", eps=.0001), lambda=0)
mcp <- coef(ncvsurv(X, y, lambda.min=0,penalty="MCP",eps=.0001), lambda=0)
expect_equivalent(scad, beta, tolerance=.01)
expect_equivalent(mcp,beta, tolerance=.01)


# ncvsurv works for cox regression ----------------------------------------

n <- 100
p <- 10
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
beta <- coef(coxph(y ~ X))
scad <- coef(ncvsurv(X, y, lambda.min=0, penalty="SCAD", eps=.0001), lambda=0)
mcp <- coef(ncvsurv(X, y, lambda.min=0,penalty="MCP",eps=.0001), lambda=0)
expect_equivalent(scad, beta, tolerance=.01)
expect_equivalent(mcp, beta, tolerance=.01)


# ncvsurv agrees with coxnet ----------------------------------------------

library(glmnet)
n <- 100
p <- 25
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
par(mfrow=c(2,1))

nlasso <- coef(nfit <- ncvsurv(X, y, penalty="lasso", lambda.min=0.001))
plot(nfit, log=TRUE)
glasso <- as.matrix(coef(gfit <- glmnet(X, y, family="cox", lambda=nfit$lambda)))
plot(gfit, "lambda")
expect_equivalent(nlasso, glasso, tolerance=.02)

expect_equivalent(predict(nfit, X, "link"), predict(gfit, X, type="link"), tolerance=.01)
expect_equivalent(predict(nfit, X, "response"), predict(gfit, X, type="response"), tolerance=.05)


# lasso/scad/mcp all seem to work: ncvsurv --------------------------------

n <- 100
p <- 10
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
par(mfrow=c(3,1))

fit <- ncvsurv(X, y, penalty="lasso")
plot(fit, main="lasso")
fit <- ncvsurv(X, y, penalty="SCAD")
plot(fit, main="scad")
fit <- ncvsurv(X, y, penalty="MCP")
plot(fit, main="mcp")


# logLik() is correct -----------------------------------------------------

n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
y <- Surv(rexp(n), rbinom(n, 1, 0.5))

fit.mle <- coxph(y~X)
fit <- ncvsurv(X, y, lambda.min=0)
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=.001)

nfit <- ncvsurv(X, y, penalty="lasso")
l <- nfit$lambda
gfit <- glmnet(X, y, family="cox", lambda=l)
-2*logLik(nfit)[50]
glmnet:::coxnet.deviance(predict(gfit, X, s=l[50]), y)


# penalty.factor seems to work: ncvsurv -----------------------------------

n <- 50
p <- 4
X <- matrix(rnorm(n*p), ncol=p)
X <- prcomp(X)$x
y <- Surv(rexp(n, 1), rbinom(n, 1, 0.5))
penalty.factor=c(0,0,1,10)

par(mfrow=c(2,1))
fit <- ncvsurv(X, y)
plot(fit)
fit <- ncvsurv(X, y, penalty.factor=penalty.factor)
plot(fit)


# ncvsurv dependencies work -----------------------------------------------

# Predict
fit <- ncvsurv(X, y)
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')
S <- predict(fit, X, 'survival', lambda=0.04)
plot(S)
S <- predict(fit, X[1,], 'survival', lambda=0.04)
plot(S)
p <- predict(fit, X, 'median')

# Coersion
fit <- ncvsurv(as.data.frame(X), y)

# Loss
eta <- predict(fit, X, 'link', lambda=0.1)
ncvreg:::loss.ncvsurv(y, eta)

# Summary
summary(fit, which=10)

# Linear predictors + residuals
fit <- ncvsurv(X, y, lambda.min=0, eps=1e-12)
sfit <- coxph(y~X)
expect_equivalent(fit$linear.predictors[, 100], sfit$linear.predictors[order(y)])
r1 <- residuals(fit, lambda=0)
r2 <- residuals(sfit, type='deviance')
expect_equivalent(residuals(fit, lambda=0), residuals(sfit, type='deviance'), tolerance=0.1)


# cox survival predictions agree with Kaplan-Meier ------------------------

n <- 50
p <- 5
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, rep(0, p-2))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.75))

fit <- ncvsurv(X, y)
S <- predict(fit, X[1,], which=1, type='survival')
km <- survfit(y~1)
plot(km, conf.int=FALSE, mark.time=FALSE, xlim=c(0,10), lwd=10, col="gray")
lines(fit$time, S(fit$time), type="s", col="slateblue", lwd=2)
median(km)
predict(fit, X[1,], which=1, type='median')


# cox survival predictions agree with coxph -------------------------------

n <- 50
p <- 5
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, rep(0, p-2))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
df <- data.frame(y, X)

fit <- ncvsurv(X, y, lambda.min=0)
sfit <- coxph(y~., data=df)
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (i in 1:9) {
  S <- predict(fit, X[i,], which=100, type='survival')
  km.cox <- survfit(sfit, newdata = df[i,], type="kalbfleisch-prentice")
  plot(km.cox, conf.int=FALSE, mark.time=FALSE, xlim=c(0,10), lwd=10, col="gray")
  lines(fit$time, S(fit$time), type="s", col="slateblue", lwd=2)
}


# cox survival predictions work -------------------------------------------

data(Lung, package='ncvreg')
X <- Lung$X
y <- Lung$y

fit <- ncvsurv(X, y)
M <- predict(fit, X, type='median')
M[1:10, 1:10]
S <- predict(fit, X, lambda=0.3, type='survival')
sapply(S, function(f) f(100))
plot(S, xlim=c(0,200))
S <- predict(fit, X, lambda=0.05, type='survival')
plot(S, xlim=c(0,200))
S <- predict(fit, X, lambda=0.44, type='survival')
plot(S, xlim=c(0,200))


# Error handling ----------------------------------------------------------

data(Lung, package='ncvreg')
X <- as.matrix(Lung$X)
storage.mode(X) <- 'integer'
y <- Lung$y
expect_silent(ncvsurv(X, y, nlambda=2, lambda.max=0.99))
expect_error(ncvsurv('a', rnorm(10)))
expect_error(ncvsurv(X, y, gamma=0.5))
expect_error(ncvsurv(X, y, gamma=1.5, penalty='SCAD'))
expect_error(ncvsurv(X, y, nlambda=1))
expect_error(ncvsurv(X, y, alpha=0))
expect_error(ncvsurv(X, y, penalty.factor=c(1,2)))
expect_error(ncvsurv(X, c(NA, y[-1])))
