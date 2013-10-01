predict.cv.ncvreg <- function(object, X, lambda=object$lambda.min, which, type=c("link","response","class","coefficients","vars","groups","norm"), ...) {
  type <- match.arg(type)
  predict.ncvreg(object$fit, X=X, lambda=lambda, which=which, type=type, ...)
}
coef.cv.ncvreg <- function(object, which=object$min, ...) {
  coef.ncvreg(object$fit, which=which, ...)
}
