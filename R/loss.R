loss.ncvreg <- function(y, yhat, family) {
  n <- length(y)
  if (family=="gaussian") {
    val <- (y-yhat)^2
  } else if (family=="binomial") {
    if (is.matrix(yhat)) {
      val <- matrix(NA, nrow=nrow(yhat), ncol=ncol(yhat))
      val[y==1,] <- -2*log(yhat[y==1, , drop=FALSE])
      val[y==0,] <- -2*log(1-yhat[y==0, , drop=FALSE])      
    } else {
      val <- numeric(length(y))
      val[y==1] <- -2*log(yhat[y==1])
      val[y==0] <- -2*log(1-yhat[y==0])
    }
  } else if (family=="poisson") {
    yly <- y*log(y)
    yly[y==0] <- 0
    val <- 2*(yly - y + yhat - y*log(yhat))
  }
  
  val
}
