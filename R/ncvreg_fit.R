ncvreg_fit <- function(X, y, family=c("gaussian","binomial"), penalty=c("MCP", "SCAD", "lasso"), gamma=3, alpha=1, lambda, eps=.001, max.iter=1000, dfmax=p+1, penalty.factor=rep(1, ncol(X)), warn=TRUE) {
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  n <- length(y)
  p <- ncol(X)
  nlam <- length(lambda)
  if (family=="gaussian") {
    fit <- .C("cdfit_gaussian", double(p*nlam), double(nlam), integer(nlam), as.double(X), as.double(y), as.integer(n), as.integer(p), penalty, as.double(lambda), as.integer(nlam), as.double(eps), as.integer(max.iter), as.double(gamma), as.double(penalty.factor), as.double(alpha), as.integer(dfmax), as.integer(TRUE))
    b <- matrix(fit[[1]], nrow=p)
    loss <- fit[[2]]
    iter <- fit[[3]]
  } else if (family=="binomial") {
    fit <- .C("cdfit_binomial", double(nlam), double(p*nlam), double(nlam), integer(nlam), as.double(X), as.double(y), as.integer(n), as.integer(p), penalty, as.double(lambda), as.integer(nlam), as.double(eps), as.integer(max.iter), as.double(gamma), as.double(penalty.factor), as.double(alpha), as.integer(dfmax), as.integer(TRUE), as.integer(warn))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    loss <- fit[[3]]
    iter <- fit[[4]]
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
