AUC.cv.ncvsurv <- function(obj, ...) {
  if (!("Y" %in% names(obj))) stop("Must run cv.ncvsurv with 'returnY=TRUE' in order to calculate AUC")
  X <- obj$Y
  y <- cbind(obj$fit$time, obj$fit$fail)
  X2 <- cbind(y, X)
  X2 <- X2[order(X2[,1]),]
  rt <- nrow(X2):1
  idx <- which(X2[,2] == 1)
  nf <- length(idx)
  n <- ncol(X)
  auct <- matrix(NA, nrow = nf, ncol = n)
  for (i in 1:nf){
    for (j in 1:n){
      xd <- X2[,1][idx[i]]
      xa <- X2[,j+2][idx[i]]
      u <- X2[,j+2][X2[,1] > xd]
      nu <- sum(u < xa)
      v <- X2[,j+2][X2[,1] > xd]
      nv <- 0.5*sum(v == xa)
      a <- (nu + nv)/(rt[idx[i]] - 1)
      auct[i,j] <- a
    }
  }
  auct[is.nan(auct)] <- 0
  as.vector(crossprod(rt[idx]-1, auct)/sum(rt[idx]-1))
}
AUC <- function(obj, ...) UseMethod("AUC")
