#' Extract Log-Likelihood
#'
#' Extract the log-likelihood of an `ncvreg` or `ncvsurv` object.
#' 
#' @param object   An `ncvreg` or `ncvsurv` object, as obtained from `ncvreg()` or `ncvsurv()`
#' @param REML     As in `logLik.lm()`
#' @param ...      For S3 compatibility
#'
#' @seealso `logLik()`
#' 
#' @rdname logLik.ncvreg
#' @export
logLik.ncvreg <- function(object, REML=FALSE, ...) {
  n <- as.double(object$n)
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
  structure(l, df=df, nobs=n, class="logLik")
}

#' @rdname logLik.ncvreg
#' @export
logLik.ncvsurv <- function(object, ...) {
  n <- as.double(object$n)
  df <- predict(object, type="nvars")
  structure(-1*object$loss, df=df, nobs=n, class="logLik")
}
