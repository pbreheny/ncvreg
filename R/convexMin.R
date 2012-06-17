convexMin <- function(beta,X,penalty,gamma,l2,family)
  {
    n <- nrow(X)
    p <- ncol(X)
    l <- ncol(beta)
    
    if (penalty=="MCP") k <- 1/gamma
    else if (penalty=="SCAD") k <- 1/(gamma-1)

    val <- NULL
    for (i in 1:l)
      {
        if (i==1) A1 <- rep(1,p)
        else A1 <- beta[-1,i]==0
        if (i==l)
          {
            L2 <- l2[i]
            U <- A1
          }
        else
          {
            if (is.na(beta[1,i+1])) break
            A2 <- beta[-1,i+1]==0
            U <- A1&A2
            L2 <- l2[i+1]
          }
        if (sum(!U)==0) next
        Xu <- X[,!U]
        if (family=="gaussian")
          {
            if (any(A1!=A2))
              {
                cmin <- min(eigen(crossprod(Xu)/n)$values)
              }
            eigen.min <- cmin - k + L2
          }
        if (family=="binomial")
          {
            if (i==l) eta <- beta[1,i] + X%*%beta[-1,i]
            else eta <- beta[1,i+1] + X%*%beta[-1,i+1]
            pi. <- exp(eta)/(1+exp(eta))
            w <- as.numeric(pi.*(1-pi.))
            w[eta > log(.9999/.0001)] <- .0001
            w[eta < log(.0001/.9999)] <- .0001
            Xu <- sqrt(w) * cbind(1,Xu)
            xwxn <- crossprod(Xu)/n
            eigen.min <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*(k-L2))))$values)
          }
        if (eigen.min < 0)
          {
            val <- i
            break
          }
      }
    return(val)
  }
