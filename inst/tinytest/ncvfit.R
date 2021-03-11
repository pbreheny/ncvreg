library(glmnet, quietly=TRUE)

# Agrees with linear regression
n <- 50
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- rnorm(p)
y <- rnorm(n, X%*%b)
fit <- lm(y~X)
m <- rep(0:1, c(1, p))
nfit <- ncvfit(cbind(1, X), y, penalty.factor=m, lambda=0)
expect_equivalent(coef(fit), nfit$beta, tol=0.001)

# Agrees with ncvreg for std matrix
n <- 50
p <- 10
L <- 30
X <- std(matrix(rnorm(n*p), n, p))
b <- rnorm(p)
y <- rnorm(n, X%*%b)
fit <- ncvreg(X, y)
m <- rep(0:1, c(1, p))
nfit <- ncvfit(cbind(1, X), y, penalty.factor=m, lambda=fit$lambda[L+1], init=coef(fit, which=L))
expect_equivalent(coef(fit, which=L+1), nfit$beta, tolerance=.001)
expect_equivalent(y - predict(fit, X, which=L+1), nfit$resid, tolerance=.001)
r <- y - predict(fit, X, which=L)
nfit <- ncvfit(cbind(1, X), y, penalty.factor=m, lambda=fit$lambda[L+1], init=coef(fit, which=L), r=r)
expect_equivalent(coef(fit, which=L+1), nfit$beta, tolerance=.001)
expect_equivalent(y - predict(fit, X, which=L+1), nfit$resid, tolerance=.001)
nfit <- ncvfit(cbind(1, X), y, penalty.factor=m, lambda=fit$lambda[L+1], init=coef(fit, which=L), r=r, xtx=c(n, rep(1, p)))
expect_equivalent(coef(fit, which=L+1), nfit$beta, tolerance=.001)
expect_equivalent(y - predict(fit, X, which=L+1), nfit$resid, tolerance=.001)

# Agrees with glmnet for nonstd matrix
n <- 50
p <- 10
L <- 30
X <- matrix(rnorm(n*p), n, p)
b <- rnorm(p)
y <- rnorm(n, X%*%b)
fit <- glmnet(X, y, standardize=FALSE)
m <- rep(0:1, c(1, p))
nfit <- ncvfit(cbind(1, X), y, penalty='lasso', penalty.factor=m, lambda=fit$lambda[L+1], init=coef(fit)[,L])
expect_equivalent(coef(fit)[,L+1], nfit$beta, tolerance=.001)
expect_equivalent(drop(y - predict(fit, X, s=fit$lambda[L+1])), nfit$resid, tolerance=.001)
r <- drop(y - predict(fit, X, s=fit$lambda[L]))
nfit <- ncvfit(cbind(1, X), y, penalty='lasso', penalty.factor=m, lambda=fit$lambda[L+1], init=coef(fit)[,L], r=r)
expect_equivalent(coef(fit)[,L+1], nfit$beta, tolerance=.001)
expect_equivalent(drop(y - predict(fit, X, s=fit$lambda[L+1])), nfit$resid, tolerance=.001)
nfit <- ncvfit(cbind(1, X), y, penalty='lasso', penalty.factor=m, lambda=fit$lambda[L+1], init=coef(fit)[,L], r=r, xtx=c(n, apply(X, 2, crossprod)/n))
expect_equivalent(coef(fit)[,L+1], nfit$beta, tolerance=.001)
expect_equivalent(drop(y - predict(fit, X, s=fit$lambda[L+1])), nfit$resid, tolerance=.001)

# Agrees with closed form in noiseless case
X <- diag(10)
y <- 1:10
expect_equivalent(ncvfit(X, y, lambda=0.2, penalty='lasso')$beta, c(0, 0, 1:8))
