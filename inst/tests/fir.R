set.seed(1)

# perm.ncvreg -------------------------------------------------------------

# Permutation of outcome
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)
pmfit <- perm.ncvreg(X, y)
pmfit <- perm.ncvreg(X, y>0, family='binomial')
pmfit <- perm.ncvreg(X, rank(y), family='poisson')

# Permutation of residuals
pmfit <- perm.ncvreg(X, y, permute='residuals')


# fir (analytic) ----------------------------------------------------------

# Linear
fit <- ncvreg(X, y)
fir(fit)
op <- par(mfrow=2:1)
plot(fir(fit))
plot(fir(fit), type="EF")
par(op)

# Logistic
n <- 50
p <- 10
X <- matrix(rnorm(n*p), n, p)
y <- rbinom(n, 1, 0.5)
fit <- ncvreg(X, y, lambda.min=0, family='binomial')
fir(fit)
op <- par(mfrow=2:1)
plot(fir(fit))
plot(fir(fit), type="EF")
par(op)

# Cox
X <- matrix(rnorm(50*10), 50, 10)
y <- cbind(rexp(50, exp(X[,1])), sample(rep(0:1, c(10,40))))
fit <- ncvsurv(X, y, lambda.min=0)
fir(fit)
op <- par(mfrow=2:1)
plot(fir(fit))
plot(fir(fit), type="EF")
par(op)

# Linear w/ penalty.factor
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50, X[,10])
fit <- ncvreg(X, y, penalty.factor=rep(0:1, c(3, 7)))
fir(fit)
op <- par(mfrow=2:1)
plot(fir(fit))
plot(fir(fit), type="EF")
par(op)
