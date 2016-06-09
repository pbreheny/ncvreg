set.seed(1)
n <- 20
p <- 5
l <- 5

.test = "std standardizes correctly"
X <- matrix(rnorm(n*p),ncol=p)
XX <- std(X)
check(apply(XX,2,mean), rep(0,5))
check(apply(XX,2,crossprod), rep(20,5))
