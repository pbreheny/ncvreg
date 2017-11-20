###################################################
.test = "ncvreg_raw works for linear regression" ##
###################################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- rnorm(p)
y <- rnorm(n, X%*%b)
# with intercept
beta <- lm(y~X)$coef
scad <- coef(ncvreg_raw(X,y,lambda=0,penalty="SCAD",eps=.0001))
mcp <- coef(ncvreg_raw(X,y,lambda=0,penalty="MCP",eps=.0001))
lasso <- coef(ncvreg_raw(X,y,lambda=0,penalty="lasso",eps=.0001))
check(scad,beta,tolerance=.01,check.attributes=FALSE)
check(mcp,beta,tolerance=.01,check.attributes=FALSE)
check(lasso,beta,tolerance=.01,check.attributes=FALSE)

# without intercept
beta0 <- lm(y~X+0)$coef
scad0 <- coef(ncvreg_raw(X,y,lambda=0,penalty="SCAD",eps=.0001,intercept=FALSE))
mcp0 <- coef(ncvreg_raw(X,y,lambda=0,penalty="MCP",eps=.0001,intercept=FALSE))
lasso0 <- coef(ncvreg_raw(X,y,lambda=0,penalty="lasso",eps=.0001,intercept=FALSE))
check(scad0,beta0,tolerance=.01,check.attributes=FALSE)
check(mcp0,beta0,tolerance=.01,check.attributes=FALSE)
check(lasso0,beta0,tolerance=.01,check.attributes=FALSE)

###############################################
.test = "logLik() is correct for ncvreg_raw" ##
###############################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)

# with intercept
y <- rnorm(n)
fit.mle <- lm(y~X)
fit <- ncvreg_raw(X, y, lambda.min=0)
check(logLik(fit)[100],  logLik(fit.mle)[1], check.attributes=FALSE, tol= .001)
check(AIC(fit)[100],  AIC(fit.mle), check.attributes=FALSE, tol= .001)

# without intercept
y <- rnorm(n)
fit0.mle <- lm(y~X+0)
fit0 <- ncvreg_raw(X, y, intercept = FALSE, lambda.min=0)
check(logLik(fit0)[100],  logLik(fit0.mle)[1], check.attributes=FALSE, tol= .001)
check(AIC(fit0)[100],  AIC(fit0.mle), check.attributes=FALSE, tol= .001)

################################################################
.test = "cv.ncvreg() works for raw data fitted by ncvreg_raw" ##
################################################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=1)

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y, raw = TRUE)
plot(cvfit, type="all")
print(summary(cvfit))
print(predict(cvfit, type="coefficients"))
print(predict(cvfit, type="vars"))
print(predict(cvfit, type="nvars"))

b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.ncvreg(X, y, raw = TRUE)
plot(cvfit, type="all")

b <- rep(0, 10)
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.ncvreg(X, y, raw = TRUE)
plot(cvfit, type="all")

#########################################
.test = "ncvreg_raw dependencies work" ##
#########################################

# Predict
fit <- ncvreg_raw(X, y, lambda.min=0)
p <- predict(fit, X, 'link', lambda=0.1)
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')

# Integers
X <- matrix(rpois(500, 1), 50, 10)
y <- rpois(50, 1)
fit <- ncvreg_raw(X, y)

# Data frame
fit <- ncvreg_raw(as.data.frame(X), y)

# Penalty factor
fit <- ncvreg_raw(as.data.frame(X), y, penalty.factor=c(0:9))

# User lambdas
fit <- ncvreg_raw(as.data.frame(X), y, lambda=c(1, 0.1, 0.01))

# ReturnX
fit <- ncvreg_raw(as.data.frame(X), y, returnX=TRUE)

# Constant columns
fit <- ncvreg_raw(cbind(5, X), y)

# Plot
plot(fit)
plot(fit, log.l=TRUE)

# Summary
summary(fit, which=10)
summary(fit, lam=0.05)
