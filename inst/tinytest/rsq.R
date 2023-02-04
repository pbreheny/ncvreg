library(survival, quietly=TRUE)

# Linear
data(Prostate, package='ncvreg')
X <- Prostate$X
y <- Prostate$y
cvfit <- cv.ncvreg(X, y)
summary(cvfit)
plot(cvfit, type='rsq')
summary(lm(y~X))

# Residuals: Gaussian
lmfit <- lm(y~X)
fit <- ncvreg(X, y, lambda.min=0)
nl <- c(nrow(X), length(fit$lambda))
expect_equal(dim(fit$linear.predictors), nl)
expect_equal(dim(residuals(fit)), nl)
expect_equal(fit$linear.predictors[,100], lmfit$fitted.values)
expect_equal(residuals(fit)[,100], residuals(lmfit))

# Logistic
data(Heart, package='ncvreg')
X <- Heart$X
y <- Heart$y
cvfit <- cv.ncvreg(X, y, family='binomial')
summary(cvfit)
plot(cvfit, type='rsq')
l1 = logLik(glm(y~X, family='binomial'))[1]
l0 = logLik(glm(y~1, family='binomial'))[1]
1 - exp(-2/length(y) * (l1 - l0))

# Residuals: Logistic
glmfit <- glm(y~X, family='binomial')
fit <- ncvreg(X, y, family='binomial', lambda.min=0)
nl <- c(nrow(X), length(fit$lambda))
expect_equal(dim(fit$linear.predictors), nl)
expect_equal(dim(residuals(fit)), nl)
expect_equal(fit$linear.predictors[,100], glmfit$linear.predictors, tolerance=0.0001)
expect_equal(residuals(fit)[,100], residuals(glmfit), tolerance=0.0001)

# Residuals: Poisson
y <- rpois(nrow(X), 1)
glmfit <- glm(y~X, family='poisson')
fit <- ncvreg(X, y, family='poisson', lambda.min=0)
nl <- c(nrow(X), length(fit$lambda))
expect_equal(dim(fit$linear.predictors), nl)
expect_equal(dim(residuals(fit)), nl)
expect_equal(fit$linear.predictors[,100], glmfit$linear.predictors, tolerance=0.0001)
expect_equal(residuals(fit)[,100], residuals(glmfit), tolerance=0.0001)

# Cox
data(Lung, package='ncvreg')
X <- Lung$X
y <- Lung$y
cvfit <- cv.ncvsurv(X, y)
summary(cvfit)
plot(cvfit, type='rsq')
summary(coxph(y~X))

# Cox residuals are checked in ncvsurv.R
