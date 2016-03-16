#source("~/dev/.ncvreg.setup.R")
require(parallel)
cl <- makeCluster(4)
#clusterCall(cl, function() source("~/dev/.ncvreg.setup.R"))

################################################
.test = "parallel version of cv.ncvreg works" ##
################################################
n <- 40
p <- 100
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- rnorm(n, mean=X%*%b, sd=2)
yb <- y > .5
yp <- rpois(n, exp(X%*%b/3))
par(mfrow=c(2,2))

cvfit <- cv.ncvreg(X, y, penalty="lasso", nfolds=2, seed=1)
plot(cvfit)
cvfit <- cv.ncvreg(X, y, penalty="lasso", nfolds=2, seed=1)
plot(cvfit)
cvfit <- cv.ncvreg(X, y, penalty="lasso", nfolds=2, seed=1, cluster=cl)
plot(cvfit)
cvfit <- cv.ncvreg(X, y, penalty="lasso", nfolds=2, seed=1, cluster=cl)
plot(cvfit)

n <- 40
p <- 5000
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- rnorm(n, mean=X%*%b, sd=2)
print(system.time(fit <- cv.ncvreg(X, y, penalty="lasso", nfolds=n, seed=1)))
print(system.time(fit <- cv.ncvreg(X, y, penalty="lasso", nfolds=n, seed=1, cluster=cl)))
stopCluster(cl)

#################################################
.test = "parallel version of cv.ncvsurv works" ##
#################################################
cl <- makeCluster(4)
n <- 50
p <- 100
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- cbind(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
par(mfrow=c(2,2))

cvfit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=2, seed=1)
plot(cvfit)
cvfit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=2, seed=2)
plot(cvfit)
cvfit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=2, seed=1, cluster=cl)
plot(cvfit)
cvfit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=2, seed=2, cluster=cl)
plot(cvfit)

n <- 40
p <- 5000
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- cbind(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
print(system.time(fit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=n, seed=1)))
print(system.time(fit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=n, seed=1, cluster=cl)))
stopCluster(cl)
