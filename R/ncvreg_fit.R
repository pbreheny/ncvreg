ncvreg_fit <- function(X, y, family=c("gaussian","binomial","poisson"), penalty=c("MCP", "SCAD", "lasso"), gamma=3, alpha=1, lambda, eps=.001, max.iter=1000, dfmax=p+1, penalty.factor=rep(1, ncol(X)), warn=TRUE) {
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  n <- length(y)
  p <- ncol(X)
  nlam <- length(lambda)
  
  if (family=="gaussian") {
    res <- .Call("cdfit_gaussian", X, as.numeric(y), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), TRUE)
    b <- matrix(res[[1]], p, nlam)
    loss <- res[[2]]
    iter <- res[[3]]
  } else if (family=="binomial") {
    res <- .Call("cdfit_binomial", X, as.numeric(y), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), TRUE, as.integer(warn))
    b0 <- numeric(nlam)
    b <- rbind(res[[1]], matrix(res[[2]], p, nlam))
    loss <- res[[3]]
    iter <- res[[4]]
  } else if (family=="poisson") {
    res <- .Call("cdfit_poisson", X, as.numeric(y), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), TRUE, as.integer(warn))
    b0 <- numeric(nlam)
    b <- rbind(res[[1]], matrix(res[[2]], p, nlam))
    loss <- res[[3]]
    iter <- res[[4]]
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
