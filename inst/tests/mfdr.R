set.seed(1)

# Linear ------------------------------------------------------------------

# Low dimensional
n <- 50
p <- 10
X <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
fit <- ncvreg(X, y)
head(mfdr(fit))

op <- par(mfrow=2:1)
plot(mfdr(fit))
plot(mfdr(fit), type="EF")
par(op)

# High dimensional
n <- 50
p <- 100
X <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
fit <- ncvreg(X, y)
head(mfdr(fit))

op <- par(mfrow=2:1)
plot(mfdr(fit))
plot(mfdr(fit), type="EF")
par(op)


# Logistic ----------------------------------------------------------------

# Low dimensional
n <- 50
p <- 10
X <- matrix(rnorm(n*p), n, p)
y <- rbinom(n, 1, binomial()$linkinv(apply(X[,1:4], 1, sum)))
fit <- ncvreg(X, y, lambda.min=0, family='binomial', returnX=TRUE)
head(mfdr(fit), n=10)
summary(fit, which=10)

op <- par(mfrow=2:1)
plot(mfdr(fit))
plot(mfdr(fit), type="EF")
par(op)

# High dimensional
n <- 50
p <- 100
X <- matrix(rnorm(n*p), n, p)
y <- rbinom(n, 1, 0.5)
fit <- ncvreg(X, y, family='binomial', returnX=TRUE)
head(mfdr(fit), n=10)

op <- par(mfrow=2:1)
plot(mfdr(fit))
plot(mfdr(fit), type="EF")
par(op)


# Cox ---------------------------------------------------------------------

# Low dimensional
X <- matrix(rnorm(50*10), 50, 10)
y <- cbind(rexp(50, exp(X[,1])), sample(rep(0:1, c(10,40))))
fit <- ncvsurv(X, y, lambda.min=0, returnX=TRUE)
head(mfdr(fit), n=10)
summary(fit, which=10)

op <- par(mfrow=2:1)
plot(mfdr(fit))
plot(mfdr(fit), type="EF")
par(op)

#############################################################
.test = "mfdr works for Cox regression when X is supplied" ##
#############################################################
m1 <- mfdr(fit)
fit <- ncvsurv(X, y, lambda.min=0)
m2 <- mfdr(fit, X)
check(m1$EF, m2$EF)


# Cox: HD
X <- matrix(rnorm(100*100), 100, 100)
y <- cbind(rexp(100, exp(X[,1])), sample(rep(0:1, c(10,40))))
fit <- ncvsurv(X, y, returnX=TRUE)
mfdr(fit)

op <- par(mfrow=2:1)
plot(mfdr(fit))
plot(mfdr(fit), type="EF")
par(op)


# Penalty factors ---------------------------------------------------------

# Linear
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50, X[,10])
fit <- ncvreg(X, y, penalty.factor=rep(0:1, c(3, 7)))
mfdr(fit)

op <- par(mfrow=2:1)
plot(mfdr(fit))
plot(mfdr(fit), type="EF")
par(op)

# Local -------------------------------------------------------------------

# High dimensional
n <- 50
p <- 100
X <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
fit <- ncvreg(X, y)
local_mfdr(fit, 0.1, method="ashr")
local_mfdr(fit, 0.1, method="kernel")
