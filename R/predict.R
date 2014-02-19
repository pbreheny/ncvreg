predict.ncvreg <- function(object, X, lambda, which=1:length(object$lambda), type=c("link", "response", "class", "coefficients", "vars", "nvars"),...) {
  type <- match.arg(type)
  beta <- coef.ncvreg(object, lambda=lambda, which=which, drop=FALSE)
  if (type=="coefficients") return(beta)
  if (type=="nvars") return(apply(beta[-1,,drop=FALSE]!=0,2,sum))
  if (type=="vars") return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which)))
  eta <- sweep(X %*% beta[-1,,drop=FALSE], 2, beta[1,], "+")
  if (object$family=="gaussian" | type=="link") return(drop(eta))
  mu <- if (object$family=="binomial") exp(eta)/(1+exp(eta)) else if (object$family=="poisson") exp(eta)
  if (type=="response") return(drop(mu))
  if (type=="class") {
    if (object$family=="binomial") {
      return(drop(1*(eta>0)))
    } else {
      stop("type='class' can only be used with family='binomial'")
    }
  }
}
coef.ncvreg <- function(object, lambda, which=1:length(object$lambda), drop=TRUE, ...) {
  if (!missing(lambda)) {
    ind <- approx(object$lambda,seq(object$lambda),lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    beta <- (1-w)*object$beta[,l,drop=FALSE] + w*object$beta[,r,drop=FALSE]
    if (length(lambda) > 1) colnames(beta) <- round(lambda,4)
  }
  else beta <- object$beta[,which,drop=FALSE]
  if (drop) return(drop(beta)) else return(beta)
}
