setupLambda <- function(X, y, family, alpha, lambda.min, nlambda, penalty.factor)
{
  n <- nrow(X)
  p <- ncol(X)
  
  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    fit <- glm(y~X[, -ind], family=family)
  } else fit <- glm(y~1, family=family)
  if (family=="gaussian") {
    z <- crossprod(X[,ind, drop=FALSE], fit$residuals) / n
  } else {
    z <- crossprod(X[,ind, drop=FALSE], residuals(fit, "working") * fit$weights) / n
  }
  lambda.max <- max(abs(z)/(alpha*penalty.factor[ind]))
  
  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  lambda
}
