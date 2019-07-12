predict.ncvsurv <- function(object, X, type=c("link", "response", "survival",
                                              "median", "coefficients", "vars",
                                              "nvars"),
                            lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)
  if (type %in% c("coefficients", "vars", "nvars")) {
    return(predict.ncvreg(object=object, X=X, type=type, lambda=lambda, which=which, ...))
  }
  if (!missing(lambda)) {
    ind <- approx(object$lambda,seq(object$lambda),lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    x <- ind %% 1
    beta <- (1-x)*object$beta[,l,drop=FALSE] + x*object$beta[,r,drop=FALSE]
    colnames(beta) <- lamNames(lambda)
  } else {
    beta <- object$beta[,which,drop=FALSE]
  }

  eta <- X %*% beta
  if (type=='link') return(drop(eta))
  if (type=='response') return(drop(exp(eta)))

  if (!missing(lambda)) {
    W <- (1-x)*exp(object$Eta)[,l,drop=FALSE] + x*exp(object$Eta)[,r,drop=FALSE]
  } else {
    W <- exp(object$Eta)[,which,drop=FALSE]
  }
  if (type == 'survival' & ncol(W) > 1) stop('Can only return type="survival" for a single lambda value')
  if (type == 'survival') val <- vector('list', length(eta))
  if (type == 'median') val <- matrix(NA, nrow(eta), ncol(eta))
  for (j in 1:ncol(eta)) {
    # Estimate baseline hazard
    w <- W[,j]
    r <- rev(cumsum(rev(w)))
    a <- ifelse(object$fail, (1-w/r)^(1/w), 1)
    S0 <- c(1, cumprod(a))
    x <- c(0, object$time)
    for (i in 1:nrow(eta)) {
      S <- S0^exp(eta[i,j])
      if (type == 'survival') val[[i]] <- approxfun(x, S, method='constant', ties=function(x) tail(x, 1))
      if (type == 'median') {
        if (any(S < 0.5)) {
          val[i,j] <- x[min(which(S < .5))]
        }
      }
    }
  }
  if (type == 'survival') {
    if (nrow(eta)==1) val <- val[[1]]
    class(val) <- c('ncvsurv.func', class(val))
    attr(val, 'time') <- object$time
  }
  if (type == 'median') val <- drop(val)
  val
}
