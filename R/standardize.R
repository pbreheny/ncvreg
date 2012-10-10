standardize <- function(X)
{
  n <- nrow(X)
  center <- colMeans(X)
  X.c <- sweep(X, 2, center)
  scale <- sqrt(apply(X.c,2,crossprod)/n)
  val <- sweep(X.c, 2, scale,"/")
  attr(val,"center") <- center
  attr(val,"scale") <- scale
  val
}
unstandardize <- function(b, center, scale)
{
  beta <- matrix(0, nrow=nrow(b), ncol=ncol(b))
  beta[-1,] <- b[-1,] / scale
  beta[1,] <- b[1,] - crossprod(center, beta[-1,,drop=FALSE])
  beta
}
