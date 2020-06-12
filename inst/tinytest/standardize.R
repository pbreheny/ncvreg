set.seed(1)
n <- 20
p <- 5
l <- 5

# std standardizes correctly"
X <- matrix(rnorm(n*p),ncol=p)
XX <- std(X)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

# std works on data frames"
XX <- std(airquality)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

# std works with integer matrix"
X <- matrix(as.integer(rpois(n*p, 1)), ncol=p)
XX <- std(X)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))

# std works with vector input"
x <- rnorm(n)
XX <- std(x)
expect_equivalent(apply(XX, 2, mean), rep(0, ncol(XX)))
expect_equivalent(apply(XX, 2, crossprod), rep(nrow(XX), ncol(XX)))
