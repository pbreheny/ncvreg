logLik.ncvreg <- function(object, REML=FALSE, ...) {
  n <- as.numeric(object$n)
  df <- predict(object, type="nvars") + 1
  if (object$family=="gaussian") {
    if (REML) rdf <- n-df
    else rdf <- n
    RSS <- object$loss
    l <- -n/2 * (log(2*pi) + log(RSS) - log(rdf)) - rdf/2
    df <- df + 1
  } else if (object$family=="binomial") {
    l <- -1*object$loss
  } else if (object$family=="poisson") {
    y <- object$y
    ind <- y != 0
    l <- -object$loss + sum(y[ind]*log(y[ind])) - sum(y) - sum(lfactorial(y))
  }
  
  val <- l
  attr(val,"df") <- df
  attr(val,"nobs") <- n
  class(val) <- "logLik"
  val
}
logLik.ncvsurv <- function(object, ...) {
  n <- as.numeric(object$n)
  df <- predict(object, type="nvars")
  val <- -1*object$loss
  attr(val,"df") <- df
  attr(val,"nobs") <- n
  class(val) <- "logLik"
  val
}