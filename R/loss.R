loss.ncvreg <- function(y, yhat, family)
{
  n <- length(y)
  if (family=="gaussian") val <- (y-yhat)^2
  if (family=="binomial") {
    val <- matrix(NA, nrow=nrow(yhat), ncol=ncol(yhat))
    val[y==1,] <- -2*log(yhat[y==1, , drop=FALSE])
    val[y==0,] <- -2*log(1-yhat[y==0, , drop=FALSE])
  }
  val
}
