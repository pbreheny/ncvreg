convexMin <- function(b, X, penalty, gamma, l2, family, penalty.factor, a, Delta=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  l <- ncol(b)
  
  if (penalty=="MCP") {
    k <- 1/gamma
  } else if (penalty=="SCAD") {
    k <- 1/(gamma-1)
  } else if (penalty=="lasso") {
    return(NULL)
  }
  if (l==0) return(NULL)
  
  val <- NULL
  for (i in 1:l) {
    A1 <- if (i==1) rep(1,p) else b[,i]==0
    if (i==l) {
      L2 <- l2[i]
      U <- A1
    } else {
      A2 <- b[,i+1]==0
      U <- A1&A2
      L2 <- l2[i+1]
    }
    if (sum(!U)==0) next
    Xu <- X[,!U]
    L2 <- L2*penalty.factor[!U]
    if (family=="gaussian") {
      if (any(A1!=A2)) {
        eigen.min <- min(eigen(crossprod(Xu)/n - diag(k-L2, sum(!U), sum(!U)))$values)
      }
    } else if (family=="binomial") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      pi. <- exp(eta)/(1+exp(eta))
      w <- as.numeric(pi.*(1-pi.))
      w[eta > log(.9999/.0001)] <- .0001
      w[eta < log(.0001/.9999)] <- .0001
      Xu <- sqrt(w) * cbind(1,Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*(k-L2))))$values)
    } else if (family=="poisson") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      mu <- exp(eta)
      w <- as.numeric(mu)
      Xu <- sqrt(w) * cbind(1,Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*(k-L2))))$values)
    } else if (family=="cox") {
      eta <- if (i==l) X%*%b[,i] else X%*%b[,i+1]
      haz <- exp(eta)
      rsk <- rev(cumsum(rev(haz)))
      W <- outer(drop(haz), rsk, "/")
      h <- apply(W, 1, function(x) cumsum(Delta*x*(1-x)))
      xwxn <- crossprod(sqrt(diag(h)) * Xu)/n
      eigen.min <- min(eigen(xwxn-diag(diag(xwxn)*(k-L2), nrow(xwxn), ncol(xwxn)))$values)
    }
    
    if (eigen.min < 0) {
      val <- i
      break
    }
  }
  val
}
