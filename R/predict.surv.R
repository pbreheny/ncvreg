predict.ncvsurv <- function(object, X, type=c("link", "response", "coefficients", "vars", "nvars"), 
                            lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)
  if (type!="response") return(predict.ncvreg(object=object, X=X, type=type, lambda=lambda, which=which, ...))
  beta <- coef.ncvreg(object, lambda=lambda, which=which, drop=FALSE)
  if (length(object$penalty.factor)==nrow(object$beta)) {
    eta <- X %*% beta
  } else {
    eta <- sweep(X %*% beta[-1,,drop=FALSE], 2, beta[1,], "+")
  }
  if (object$model=="cox") resp <- exp(eta)
  drop(resp)
}
