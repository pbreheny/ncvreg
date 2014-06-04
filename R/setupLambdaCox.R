setupLambdaCox <- function(X, y, Delta, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)
  
  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    res <- .Call("cdfit_cox_dh", X, y, Delta, "lasso", 0, 0.001, as.integer(100), 3, penalty.factor, 
                    alpha, as.integer(p), as.integer(TRUE), as.integer(FALSE))
  } else {
    w <- 1/(n-(1:n)+1)
    h <- cumsum(Delta*w*(1-w))
    s <- Delta - cumsum(Delta*w)
    r <- s/h
    r[h==0] <- 0
  }
  zmax <- .Call("maxprod", X, r * h, ind, penalty.factor) / n
  lambda.max <- zmax/alpha
  
  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  lambda
}
