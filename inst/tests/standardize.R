set.seed(1)
n <- 20
p <- 5
l <- 5

.test = "std standardizes correctly"
X <- matrix(rnorm(n*p),ncol=p)
XX <- std(X)
check(apply(XX, 2, mean), rep(0, ncol(XX)))
check(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

.test = "std works on data frames"
XX <- std(airquality)
check(apply(XX, 2, mean), rep(0, ncol(XX)))
check(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

.test = "std works with integer matrix"
X <- matrix(as.integer(rpois(n*p, 1)), ncol=p)
XX <- std(X)
check(apply(XX, 2, mean), rep(0, ncol(XX)))
check(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

.test = "std works with vector input"
x <- rnorm(n)
XX <- std(x)
check(apply(XX, 2, mean), rep(0, ncol(XX)))
check(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))
