ncvreg_fit <- function(X, y, family=c("gaussian","binomial"), penalty=c("MCP", "SCAD", "lasso"), gamma=3, alpha=1, lambda, eps=.001, max.iter=1000, dfmax=p+1, penalty.factor=rep(1, ncol(X)), warn=TRUE) {
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  n <- length(y)
  p <- ncol(X)
  nlam <- length(lambda)
  if (family=="gaussian") {
    b <- matrix(0, p, nlam)
    loss <- numeric(nlam)
    iter <- integer(nlam)
    .Call("cdfit_gaussian", b, loss, iter, X, as.numeric(y), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), TRUE)
  } else if (family=="binomial") {
    b0 <- numeric(nlam)
    b <- matrix(0, p, nlam)
    loss <- numeric(nlam)
    iter <- integer(nlam)
    .Call("cdfit_binomial", b0, b, loss, iter, X, as.numeric(y), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), TRUE, as.integer(warn))
    b <- rbind(b0, b)
  }
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(b[p,])
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")
  
  list(beta=b, loss=loss, iter=iter, lambda=lambda)
}
