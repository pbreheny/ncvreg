require(ncvreg)
set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001, check.attributes=FALSE)}

#### Linear regression ####

# Works
y <- cbind(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*10), 50, 10)
cvfit <- cv.ncvsurv(X, y, lambda.min=0)
plot(cvfit, type='all')
print(summary(cvfit))

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
print(summary(cvfit))

# AUC
cvfit <- cv.ncvsurv(X, y, returnY=TRUE)
AUC(cvfit)
