unstandardize <- function(b, center, scale) {
  beta <- matrix(0, nrow=nrow(b), ncol=ncol(b))
  beta[-1,] <- b[-1,] / scale
  beta[1,] <- b[1,] - crossprod(center, beta[-1,,drop=FALSE])
  beta
}
