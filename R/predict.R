predict.ncvreg <- function(object, X, type=c("link", "response", "class", "coefficients", "vars", "nvars"),
                           lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)
  beta <- coef.ncvreg(object, lambda=lambda, which=which, drop=FALSE)
  if (type=="coefficients") return(beta)
  if (class(object)[1]=='ncvreg') {
    alpha <- beta[1,]
    beta <- beta[-1,,drop=FALSE]
  } else {
    beta <- beta
  }

  if (type=="nvars") return(apply(beta!=0,2,sum))
  if (type=="vars") return(drop(apply(beta!=0, 2, FUN=which)))
  eta <- sweep(X %*% beta, 2, alpha, "+")
  if (type=="link" || object$family=="gaussian") return(drop(eta))
  resp <- switch(object$family,
                 binomial = exp(eta)/(1+exp(eta)),
                 poisson = exp(eta))
  if (type=="response") return(drop(resp))
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
    colnames(beta) <- round(lambda,4)
  }
  else beta <- object$beta[, which, drop=FALSE]
  if (drop) return(drop(beta)) else return(beta)
}
