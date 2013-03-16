loss.ncvreg <- function(y, yhat, family)
{
  n <- length(y)
  if (family=="gaussian") val <- (y-yhat)^2
  if (family=="binomial") {
    if (is.matrix(yhat)) {
      val <- matrix(NA, nrow=nrow(yhat), ncol=ncol(yhat))
      val[y==1,] <- -2*log(yhat[y==1, , drop=FALSE])
      val[y==0,] <- -2*log(1-yhat[y==0, , drop=FALSE])      
    } else {
      val <- numeric(length(y))
      val[y==1] <- -2*log(yhat[y==1])
      val[y==0] <- -2*log(1-yhat[y==0])
    }
  }
  val
}
