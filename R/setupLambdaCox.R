setupLambdaCox <- function(X, y, Delta, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)

  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    nullFit <- .Call("cdfit_cox_dh", X[, -ind, drop=FALSE], y, Delta, "lasso", 0, 0.001, as.integer(100), 3, penalty.factor[-ind],
                    alpha, as.integer(p), as.integer(TRUE), as.integer(FALSE))
    eta <- as.numeric(X[, -ind, drop=FALSE] %*% nullFit[[1]])
    rsk <- rev(cumsum(rev(exp(eta))))
    s <- Delta - exp(eta)*cumsum(Delta/rsk)
  } else {
    w <- 1/(n-(1:n)+1)
    s <- Delta - cumsum(Delta*w)
  }
  zmax <- .Call("maxprod", X, s, ind, penalty.factor) / n
  lambda.max <- zmax/alpha

  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  lambda
}
