require(ncvreg)
set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001, check.attributes=FALSE)}

#### Linear regression ####

# Works
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)
cvfit <- cv.ncvreg(X, y)
plot(cvfit, type='all')
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
plot(cvfit, type='all')
print(summary(cvfit))

#### Logistic regression ####

# Works
X <- matrix(rnorm(500), 50, 10)
y <- rbinom(50, 1, 0.5)
cvfit <- cv.ncvreg(X, y, family='binomial')
plot(cvfit, type='all')
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
plot(cvfit, type='all')
print(summary(cvfit))

#### Poisson regression ####

# Works
cvfit <- cv.ncvreg(X, y, family='poisson')
plot(cvfit, type='all')
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
plot(cvfit, type='all')
print(summary(cvfit))
