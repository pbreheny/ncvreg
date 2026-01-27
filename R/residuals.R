#' Extract residuals from a ncvreg or ncvsurv fit
#'
#' Currently, only deviance residuals are supported.
#'
#' @param object   Object of class `ncvreg` or `ncvsurv`.
#' @param lambda   Values of the regularization parameter at which residuals are requested (numeric vector). For values of lambda not in the sequence of fitted models, linear interpolation is used.
#' @param which    Index of the penalty parameter at which residuals are requested (default = all indices). If `lambda` is specified, this take precedence over `which`.
#' @param drop     By default, if a single value of lambda is supplied, a vector of residuals is returned (logical; default=`TRUE`). Set `drop=FALSE` if you wish to have the function always return a matrix (see [drop()]).
#' @param ...      Not used.
#'
#' @examples
#' data(Prostate)
#' X <- Prostate$X
#' y <- Prostate$y
#' fit <- ncvreg(X, y)
#' residuals(fit)[1:5, 1:5]
#' head(residuals(fit, lambda=0.1))
#' @export

residuals.ncvreg <- function(object, lambda, which=1:length(object$lambda), drop=TRUE, ...) {

  # Calculate matrix of residuals
  if (inherits(object, "ncvsurv")) {
    for (j in 1:length(object$lambda)) {
      h <- suppressWarnings(predict(object, which = j, type = "hazard")[[1]](object$time))
      M <- object$fail - h * exp(object$linear.predictors)
      R <- sign(M) * sqrt(-2*(M + object$fail * log(object$fail - M)))
      R[h == 0,] <- 0
      R <- R[match(1:object$n, object$order),]  # Return in original order
    }
  } else if (object$family == 'gaussian') {
    R <- object$y - object$linear.predictors
  } else if (object$family == 'binomial') {
    f <- binomial()$dev.resids
    M <- binomial()$linkinv(object$linear.predictors)
    R <- vapply(1:length(object$lambda), function(j) {dr(f, object$y, M[,j])}, double(length(object$y)))
  } else if (object$family == 'poisson') {
    f <- poisson()$dev.resids
    M <- poisson()$linkinv(object$linear.predictors)
    R <- vapply(1:length(object$lambda), function(j) {dr(f, object$y, M[,j])}, double(length(object$y)))
  } else {
    stop('Residuals not implemented for this type of ncvreg object.')
  }

  # Interpolate and return
  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    out <- (1-w)*R[, l, drop=FALSE] + w*R[, r, drop=FALSE]
    colnames(out) <- round(lambda, 4)
  } else {
    out <- R[, which, drop=FALSE]
  }
  if (drop) return(drop(out)) else return(out)
}

dr <- function(f, y, m) {
  sqrt(pmax(f(y, m, rep(1, length(y))), 0)) * ((y > m) * 2 - 1)
}
