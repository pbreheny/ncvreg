setupLambda <- function(X, y, family, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)

  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    fit <- glm(y~X[, -ind], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }
  if (family=="gaussian") {
    zmax <- .Call("maxprod", X, fit$residuals, ind, penalty.factor) / n
  } else {
    zmax <- .Call("maxprod", X, residuals(fit, "working") * fit$weights, ind, penalty.factor) / n
  }
  lambda.max <- zmax/alpha

  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  } else {
    lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  }

  if (length(ind)!=p) lambda[1] <- lambda[1] * 1.000001
  lambda
}
