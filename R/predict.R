predict.ncvreg <- function(object,X,lambda,which=1:length(object$lambda),type=c("link","response","class","coefficients"),...)
  {
    type <- match.arg(type)
    beta <- coef.ncvreg(object,lambda=lambda,which=which)
    if (type=="coefficients") return(beta)
    eta <- sweep(X %*% beta[-1,,drop=FALSE], 2, beta[1,], "+")
    if (object$family=="gaussian" | type=="link") return(eta)
    pihat <- exp(eta)/(1+exp(eta))
    if (type=="response") return(pihat)
    if (type=="class") return(eta>0)
  }
coef.ncvreg <- function(object,lambda,which=1:length(object$lambda),...)
  {
    if (!missing(lambda))
      {
        ind <- approx(object$lambda,seq(object$lambda),lambda)$y
        l <- floor(ind)
        r <- ceiling(ind)
        w <- ind %% 1
        beta <- (1-w)*object$beta[,l] + w*object$beta[,r]
        if (length(lambda) > 1) colnames(beta) <- round(lambda,4)
      }
    else beta <- object$beta[,which]
    return(beta)
  }
