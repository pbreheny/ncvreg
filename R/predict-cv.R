#' @rdname predict.ncvreg
#' @export
predict.cv.ncvreg <- function(object, X, type=c("link","response","class","coefficients","vars","nvars"), which=object$min, ...) {
  type <- match.arg(type)
  predict.ncvreg(object$fit, X=X, which=which, type=type, ...)
}


#' @rdname predict.ncvreg
#' @export
coef.cv.ncvreg <- function(object, which=object$min, ...) {
  coef.ncvreg(object$fit, which=which, ...)
}

#' @rdname predict.ncvreg
#' @export
predict.cv.ncvsurv <- function(object, X,
                               type=c("link", "response", "survival", "median", "coefficients", "vars", "nvars"),
                               which=object$min, ...) {
  type <- match.arg(type)
  predict.ncvsurv(object$fit, X=X, which=which, type=type, ...)
}
