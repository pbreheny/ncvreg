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

# Gaussian w/ unpenalized
n <- 50
p <- 20
X <- matrix(rnorm(n*p), n, p)
b <- c(1, -1, 0.5, -0.5, rep(0.25, 3), rep(-0.25, 3), rep(0, 10))
y <- rnorm(n, X%*%b)
fit <- ncvreg(X, y, returnX=TRUE, penalty.factor=c(0, 0, rep(1, 18)))
summary(fit, lambda=0.2)
fit <- ncvreg(X, y, returnX=TRUE, penalty='lasso')
summary(fit, lambda=0.2)
