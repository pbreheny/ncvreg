source("~/dev/.ncvreg.setup.R")
require(survival)

##################################################################
.test = "ncvsurv works for simple cox regression, no censoring" ##
##################################################################
n <- 10
X <- matrix(runif(n), n, 1)
y <- Surv(sort(rexp(n, X), decreasing=TRUE), rep(1,n))
X <- X - mean(X)
beta <- coef(coxph(y ~ X))
eta <- X%*%beta
w <- exp(eta)/cumsum(exp(eta))
sum(X * (1-w))
scad <- coef(ncvsurv(X, y, lambda=0, penalty="SCAD", eps=.0001))
mcp <- coef(ncvsurv(X, y, lambda=0, penalty="MCP", eps=.0001))
check(scad, beta, tolerance=.01, check.attributes=FALSE)
check(mcp, beta, tolerance=.01,check.attributes=FALSE)

####################################################
.test = "ncvsurv works for simple cox regression" ##
####################################################
n <- 30
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
y[which.max(y[,1]), 2] <- 0
y[which.min(y[,1]), 2] <- 0
X <- matrix(rnorm(n), n, 1)
beta <- coef(coxph(y ~ X))
scad <- coef(ncvsurv(X, y, lambda=0, penalty="SCAD", eps=.0001))
mcp <- coef(ncvsurv(X, y, lambda=0,penalty="MCP",eps=.0001))
check(scad, beta, tolerance=.01, check.attributes=FALSE)
check(mcp,beta, tolerance=.01, check.attributes=FALSE)

#############################################
.test = "ncvsurv works for cox regression" ##
#############################################
n <- 100
p <- 10
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
beta <- coef(coxph(y ~ X))
scad <- coef(ncvsurv(X, y, lambda=0, penalty="SCAD", eps=.0001))
mcp <- coef(ncvsurv(X, y, lambda=0,penalty="MCP",eps=.0001))
check(scad, beta, tolerance=.01, check.attributes=FALSE)
check(mcp, beta, tolerance=.01, check.attributes=FALSE)

#######################################
.test = "ncvsurv agrees with coxnet" ##
#######################################
require(glmnet)
n <- 100
p <- 10
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
par(mfrow=c(2,1))

nlasso <- coef(nfit <- ncvsurv(X, y, penalty="lasso"))
plot(nfit, log=TRUE)
glasso <- as.matrix(coef(gfit <- glmnet(X, y, family="cox", lambda=nfit$lambda)))
plot(gfit, "lambda")
check(nlasso, glasso, tolerance=.01, check.attributes=FALSE)

check(predict(nfit, X, "link"), predict(gfit, X, type="link"), tolerance=.01, check.attributes=FALSE)
check(predict(nfit, X, "response"), predict(gfit, X, type="response"), tolerance=.01, check.attributes=FALSE)

############################################
.test = "lasso/scad/mcp all seem to work" ##
############################################
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

################################
.test = "logLik() is correct" ##
################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
y <- Surv(rexp(n), rbinom(n, 1, 0.5))

fit.mle <- coxph(y~X)
fit <- ncvsurv(X, y, lambda.min=0)
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)

nfit <- ncvsurv(X, y, penalty="lasso")
l <- nfit$lambda
gfit <- glmnet(X, y, family="cox", lambda=l)
-2*logLik(nfit)[50]
glmnet:::coxnet.deviance(predict(gfit, X, s=l[50]), y)

#############################################
.test = "ncvsurv handles constant columns" ##
#############################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
X[, 3:5] <- 0
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
fit <- ncvsurv(X, y)

#######################################
.test = "cv.ncvsurv() seems to work" ##
#######################################
n <- 50
p <- 100
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
par(mfrow=c(2,2))

ncvfit <- cv.ncvsurv(X, y, penalty="lasso", lambda.min=0.2)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, penalty="lasso", events.only=TRUE, lambda.min=0.2)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, lambda.min=0.4)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, events.only=TRUE, lambda.min=0.4)
plot(ncvfit)

require(glmnet)
par(mfrow=c(2,1))
ncvfit <- cv.ncvsurv(X, y, penalty="lasso")
plot(ncvfit)
gcvfit <- cv.glmnet(X, y, family="cox", grouped=TRUE, lambda=ncvfit$lambda)
plot(gcvfit)

#########################################
.test = "penalty.factor seems to work" ##
#########################################
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

######################################
.test = "cv.ncvsurv() options work" ##
######################################
n <- 50
p <- 20
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))

par(mfrow=c(2,2))
cvfit <- cv.ncvsurv(X, y)
plot(cvfit, type="all")
par(mfrow=c(2,2))
cvfit <- cv.ncvsurv(X, y, events.only=FALSE)
plot(cvfit, type="all")
print(summary(cvfit))
print(predict(cvfit, "coefficients"))
print(predict(cvfit, "vars"))
print(predict(cvfit, "nvars"))
