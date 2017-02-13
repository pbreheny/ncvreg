require(parallel)
require(survival)
cl <- makeCluster(4)

################################################
.test = "parallel version of cv.ncvreg works" ##
################################################

# Linear
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)
cvfit <- cv.ncvreg(X, y, cluster=cl)

# Logistic
y <- rbinom(50, 1, 0.5)
cvfit <- cv.ncvreg(X, y, cluster=cl, family='binomial')


#################################################
.test = "parallel version of cv.ncvsurv works" ##
#################################################
y <- Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*10), 50, 10)
cvfit <- cv.ncvsurv(X, y, cluster=cl)

####################################
.test = "parallel speedup timing" ##
####################################
n <- 40
p <- 5000
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- rnorm(n, mean=X%*%b, sd=2)
print(system.time(fit <- cv.ncvreg(X, y, penalty="lasso", nfolds=n, seed=1)))
print(system.time(fit <- cv.ncvreg(X, y, penalty="lasso", nfolds=n, seed=1, cluster=cl)))

n <- 40
p <- 5000
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- cbind(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
print(system.time(fit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=n, seed=1)))
print(system.time(fit <- cv.ncvsurv(X, y, penalty="lasso", nfolds=n, seed=1, cluster=cl)))

stopCluster(cl)
