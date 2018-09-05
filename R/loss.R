loss.ncvreg <- function(y, yhat, family) {
  n <- length(y)
  if (family=="gaussian") {
    val <- (y-yhat)^2
  } else if (family=="binomial") {
    yhat[yhat < 0.00001] <- 0.00001
    yhat[yhat > 0.99999] <- 0.99999
    if (is.matrix(yhat)) {
      val <- matrix(NA, nrow=nrow(yhat), ncol=ncol(yhat))
      if (sum(y==1)) val[y==1,] <- -2*log(yhat[y==1, , drop=FALSE])
      if (sum(y==0)) val[y==0,] <- -2*log(1-yhat[y==0, , drop=FALSE])
    } else {
      val <- numeric(length(y))
      if (sum(y==1)) val[y==1] <- -2*log(yhat[y==1])
      if (sum(y==0)) val[y==0] <- -2*log(1-yhat[y==0])
    }
  } else if (family=="poisson") {
    yly <- y*log(y)
    yly[y==0] <- 0
    val <- 2*(yly - y + yhat - y*log(yhat))
  }
  val
}
loss.ncvsurv <- function(y, eta, total=TRUE) {
  ind <- order(y[,1])
  d <- as.numeric(y[ind,2])
  if (is.matrix(eta)) {
    eta <- eta[ind, , drop=FALSE]
    r <- apply(eta, 2, function(x) rev(cumsum(rev(exp(x)))))
  } else {
    eta <- eta[ind]
    r <- rev(cumsum(rev(exp(eta))))
  }
  if (total) {
    return(-2*(crossprod(d, eta) - crossprod(d, log(r))))
  } else {
    return(-2*(eta[d==1,] - log(r)[d==1,]))
  }
}
