# Gaussian
n <- 50
p <- 20
X <- matrix(rnorm(n*p), n, p)
b <- c(1, -1, 0.5, -0.5, rep(0.25, 3), rep(-0.25, 3), rep(0, 10))
y <- rnorm(n, X%*%b)
fit <- ncvreg(X, y, returnX=TRUE)
summary(fit, lambda=0.2)
fit <- ncvreg(X, y, returnX=TRUE, penalty='lasso')
summary(fit, lambda=0.2)

# Binomial
n <- 100
p <- 20
X <- matrix(rnorm(n*p), n, p)
b <- c(1, -1, 0.5, -0.5, rep(0.25, 3), rep(-0.25, 3), rep(0, 10))
y <- rnorm(n, X%*%b) > 0
fit <- ncvreg(X, y, family='binomial', returnX=TRUE)
summary(fit, lambda=0.075)

# Cox
n <- 100
p <- 10
y <- cbind(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
fit <- ncvsurv(X, y)
summary(fit, lambda=fit$lambda[30])

# Unpenalized + constant ------------------------------------------------------
n <- 50
p <- 19
X1 <- cbind(matrix(rnorm(n*p), n, p), 0)
b1 <- c(1, -1, 0.5, -0.5, rep(0.25, 3), rep(-0.25, 3), rep(0, 10))
f1 <- c(0, 0, rep(1, 18))
X2 <- cbind(0, matrix(rnorm(n*p), n, p))
b2 <- c(0, 1, -1, 0.5, -0.5, rep(0.25, 3), rep(-0.25, 3), rep(0, 9))
f2 <- c(1, 0, 0, rep(1, 17))

# Gaussian
y1 <- rnorm(n, X1%*%b1)
fit <- ncvreg(X1, y1, returnX=TRUE, penalty.factor=f1)
summary(fit, lambda=0.2)
y2 <- rnorm(n, X2%*%b2)
fit <- ncvreg(X2, y2, returnX=TRUE, penalty.factor=f2)
summary(fit, lambda=0.2)
fit <- ncvreg(X2, y2, returnX=FALSE, penalty.factor=f2)
summary(fit, X=X2, y=y2, lambda=0.2)

# Binomial
y1 <- rnorm(n, X1%*%b1) > 0
fit <- ncvreg(X1, y1, returnX=TRUE, penalty.factor=f1, family='binomial')
summary(fit, lambda=fit$lambda[10])
y2 <- rnorm(n, X2%*%b2) > 0
fit <- ncvreg(X2, y2, returnX=TRUE, penalty.factor=f2, family='binomial')
summary(fit, lambda=fit$lambda[10])
fit <- ncvreg(X2, y2, returnX=FALSE, penalty.factor=f2, family='binomial')
summary(fit, X=X2, y=y2, lambda=fit$lambda[10])

# Cox
y1 <- cbind(rexp(n, exp(X1%*%b1)), 1)
fit <- ncvsurv(X1, y1, penalty.factor=f1)
summary(fit, lambda=fit$lambda[10])
y2 <- cbind(rexp(n, exp(X1%*%b1)), 1)
fit <- ncvsurv(X2, y2, penalty.factor=f2)
summary(fit, lambda=fit$lambda[10])
fit <- ncvsurv(X2, y2, returnX=FALSE, penalty.factor=f2)
summary(fit, X=X2, y=y2, lambda=fit$lambda[10])
