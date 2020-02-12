setupLambdaCox <- function(X, y, Delta, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)

  # Fit to unpenalized covariates
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    nullFit <- survival::coxph(survival::Surv(y, Delta) ~ X[, -ind, drop=FALSE])
    eta <- nullFit$linear.predictors
    rsk <- rev(cumsum(rev(exp(eta))))
    s <- Delta - exp(eta)*cumsum(Delta/rsk)
  } else {
    w <- 1/(n-(1:n)+1)
    s <- Delta - cumsum(Delta*w)
  }

  # Determine lambda.max
  zmax <- .Call("maxprod", X, s, ind, penalty.factor) / n
  lambda.max <- zmax/alpha

  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
  else lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
  lambda
}
