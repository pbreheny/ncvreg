set.seed(1)
n <- 20
p <- 5
l <- 5

# std standardizes correctly
X <- matrix(rnorm(n*p),ncol=p)
XX <- std(X)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

# std works on data frames
XX <- std(airquality)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

# std works with integer matrix
X <- matrix(as.integer(rpois(n*p, 1)), ncol=p)
XX <- std(X)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

# std works with vector input
x <- rnorm(n)
XX <- std(x)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

# std works with new data
X <- matrix(rnorm(50), 10, 5)
S <- std(X)
XX <- matrix(rnorm(50), 10, 5)
SS <- std(S, X)
expect_true(all.equal(colMeans(SS), rep(0, ncol(SS)), tol=0.1))
expect_true(all.equal(colMeans(SS^2), rep(1, ncol(SS)), tol=0.1))
SS <- std(S, XX)
expect_inherits(all.equal(colMeans(SS), rep(0, ncol(SS)), tol=0.1), 'character')
expect_inherits(all.equal(colMeans(SS^2), rep(1, ncol(SS)), tol=0.1), 'character')

# Predictions all agree
y <- rnorm(10)
fit <- ncvreg(X, y)
P1 <- predict(fit, X)
fit <- ncvreg(S, y)
P2 <- predict(fit, S)
expect_equivalent(P1, P2)
y <- rnorm(10)
fit <- ncvreg(X, y)
P1 <- predict(fit, XX)
fit <- ncvreg(S, y)
P2 <- predict(fit, SS)
expect_equivalent(P1, P2)
