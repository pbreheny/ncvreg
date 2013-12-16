predict.cv.ncvreg <- function(object, type=c("link","response","class","coefficients","vars","nvars"), which=object$min, ...) {
  type <- match.arg(type)
  predict.ncvreg(object$fit, which=which, type=type, ...)
}
coef.cv.ncvreg <- function(object, which=object$min, ...) {
  coef.ncvreg(object$fit, which=which, ...)
}
