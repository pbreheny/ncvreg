require(ncvreg)
if (Sys.getenv('USER')=="pbreheny") {
  require(parallel)
  cl <- makeCluster(4)

  # Linear
  X <- matrix(rnorm(500), 50, 10)
  y <- rnorm(50)
  cvfit <- cv.ncvreg(X, y, cluster=cl)

  # Logistic
  y <- rbinom(50, 1, 0.5)
  cvfit <- cv.ncvreg(X, y, cluster=cl, family='binomial')

  # Cox
  y <- Surv(rexp(50), sample(rep(0:1, c(10,40))))
  X <- matrix(rnorm(50*10), 50, 10)
  cvfit <- cv.ncvsurv(X, y, cluster=cl)
}
