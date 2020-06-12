set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001)}

# Works
X <- matrix(rnorm(50*10), 50, 10)
y <- cbind(rexp(50, exp(X[,1])), sample(rep(0:1, c(10,40))))
cvfit <- cv.ncvsurv(X, y, lambda.min=0)
plot(cvfit)
plot(cvfit, type='all')
summary(cvfit)

# Predict
p <- predict(cvfit, X, 'link', lambda=0.1)
p <- predict(cvfit, X, 'link')
p <- predict(cvfit, X, 'response')
p <- predict(cvfit, X, 'coef')
p <- predict(cvfit, X, 'vars')
p <- predict(cvfit, X, 'nvars')
S <- predict(cvfit, X, 'survival', lambda=0.1)
plot(S)
S <- predict(cvfit, X[1,], 'survival', lambda=0.1)
plot(S)
p <- predict(cvfit, X, 'median')

# LOOCV
y <- cbind(rexp(25), sample(rep(0:1, c(5,20))))
X <- matrix(rnorm(25*5), 25, 5)
cvfit <- cv.ncvsurv(X, y, nfolds=25)
plot(cvfit, type='all')
summary(cvfit)

# AUC
cvfit <- cv.ncvsurv(X, y, returnY=TRUE)
AUC(cvfit)

############################################
# cv.ncvsurv() agrees with glmnet
############################################
require(survival)
n <- 50
p <- 100
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
par(mfrow=c(2,2))

ncvfit <- cv.ncvsurv(X, y, penalty="lasso", lambda.min=0.2)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, penalty="lasso", events.only=FALSE, lambda.min=0.2)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, lambda.min=0.4)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, events.only=FALSE, lambda.min=0.4)
plot(ncvfit)

require(glmnet)
par(mfrow=c(2,1))
gcvfit <- cv.glmnet(X, y, family="cox", grouped=TRUE, foldid=ncvfit$fold, lambda.min=0.01)
ncvfit <- cv.ncvsurv(X, y, penalty="lasso", fold=ncvfit$fold, lambda=gcvfit$lambda)
plot(ncvfit)
plot(gcvfit, sign.lambda = -1)

require(glmnet)
par(mfrow=c(2,1))
ncvfit1 <- cv.ncvsurv(X, y, penalty="lasso", fold=ncvfit$fold, se='quick')
ncvfit2 <- cv.ncvsurv(X, y, penalty="lasso", fold=ncvfit$fold, se='bootstrap')
ylim <- range(c(ncvfit1$cve+ncvfit1$cvse, ncvfit1$cve-ncvfit1$cvse, ncvfit2$cve+ncvfit2$cvse, ncvfit2$cve-ncvfit2$cvse))
plot(ncvfit1, ylim=ylim)
plot(ncvfit2, ylim=ylim)

######################################
# cv.ncvsurv() options work
######################################
n <- 50
p <- 20
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- cbind(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))

par(mfrow=c(2,2))
cvfit <- cv.ncvsurv(X, y)
plot(cvfit, type="all")
par(mfrow=c(2,2))
summary(cvfit)
predict(cvfit, type="coefficients")
predict(cvfit, type="vars")
predict(cvfit, type="nvars")

