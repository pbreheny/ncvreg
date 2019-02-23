set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001, check.attributes=FALSE)}

#### Linear regression ####

# Works
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)
cvfit <- cv.ncvreg(X, y)
print(summary(cvfit))

# Predict
b <- coef(cvfit)
p <- predict(cvfit, X, type='link')
p <- predict(cvfit, X, type='response')
p <- predict(cvfit, X, type='coef')
p <- predict(cvfit, X, type='vars')
p <- predict(cvfit, X, type='nvars')

# Integers
X <- matrix(rpois(500, 1), 50, 10)
y <- rpois(50, 1)
cvfit <- cv.ncvreg(X, y)

# Data frame
cvfit <- cv.ncvreg(as.data.frame(X), y)

# LOOCV
X <- matrix(rnorm(25*4), 25, 4)
y <- rnorm(25)
cvfit <- cv.ncvreg(X, y, nfolds=25)
print(summary(cvfit))
cvfit <- cv.ncvreg(X, y, fold=1:25)

#### Logistic regression ####

# Works
X <- matrix(rnorm(500), 50, 10)
y <- rbinom(50, 1, 0.5)
cvfit <- cv.ncvreg(X, y, family='binomial')
print(summary(cvfit))

# Predict
b <- coef(cvfit)
p <- predict(cvfit, X, type='link')
p <- predict(cvfit, X, type='response')
p <- predict(cvfit, X, type='class')
p <- predict(cvfit, X, type='coef')
p <- predict(cvfit, X, type='vars')
p <- predict(cvfit, X, type='nvars')

# LOOCV
X <- matrix(rnorm(30*2), 30, 2)
y <- rbinom(30, 1, 0.5)
cvfit <- cv.ncvreg(X, y, nfolds=30, family='binomial')
print(summary(cvfit))
cvfit <- cv.ncvreg(X, y, fold=1:30, family='binomial')

#### Poisson regression ####

# Works
cvfit <- cv.ncvreg(X, y, family='poisson')
print(summary(cvfit))

# Predict
b <- coef(cvfit)
p <- predict(cvfit, X, type='link')
p <- predict(cvfit, X, type='response')
p <- predict(cvfit, X, type='coef')
p <- predict(cvfit, X, type='vars')
p <- predict(cvfit, X, type='nvars')

# LOOCV
cvfit <- cv.ncvreg(X, y, nfolds=30, family='poisson')
print(summary(cvfit))
cvfit <- cv.ncvreg(X, y, fold=1:30, family='poisson')

######################################
.test = "cv.ncvreg() seems to work" ##
######################################
n <- 40
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- rnorm(n, mean=X%*%b, sd=2)
yb <- y > .5
yp <- rpois(n, exp(X%*%b/3))
par(mfrow=c(3,2))

require(glmnet)
gcvfit <- cv.glmnet(X, y, nfolds=n)
plot(gcvfit)
ncvfit <- cv.ncvreg(X, y, penalty="lasso", lambda=gcvfit$lambda, nfolds=n)
plot(ncvfit)
gcvfit <- cv.glmnet(X, yb, family="binomial", nfolds=n)
plot(gcvfit)
ncvfit <- cv.ncvreg(X, yb, family="binomial", penalty="lasso", lambda=gcvfit$lambda, nfolds=n)
plot(ncvfit)
cvfit <- cv.glmnet(X, yp, family="poisson")
plot(cvfit)
cvfit <- cv.ncvreg(X, yp, family="poisson", penalty="lasso", lambda=cvfit$lambda)
plot(cvfit)

##############################################
.test = "cv.ncvreg() return LP array works" ##
##############################################
n <- 100
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-3, 3, rep(0, 8))

y <- rnorm(n, mean=X%*%b, sd=1)
cvfit <- cv.ncvreg(X, y, returnY=TRUE)
cve <- apply(cvfit$Y - y, 2, crossprod)/n
check(cve, cvfit$cve, check.attributes=FALSE, tol= .001)

y <- rnorm(n, mean=X%*%b) > 0
cvfit <- cv.ncvreg(X, y, family='binomial', returnY=TRUE)
pe <- apply((cvfit$Y>=0.5)!=y, 2, mean)
check(pe, cvfit$pe, check.attributes=FALSE, tol= .001)
