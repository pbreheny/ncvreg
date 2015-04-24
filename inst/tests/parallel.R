source("~/dev/.ncvreg.setup.R")
require(parallel)
cl <- makeCluster(4)
clusterCall(cl, function() source("~/dev/.ncvreg.setup.R"))

#################################################
.test = "parallel version of cv.ncvreg workds" ##
#################################################
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
