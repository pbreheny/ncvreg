cv.ncvreg <- function(X, y, family=c("gaussian","binomial"), alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, nfolds=10, seed, trace=FALSE, ...)
  {
    ## Error checking
    family <- match.arg(family)
    if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
    if (!missing(seed)) set.seed(seed)

    ## Set up XX, yy, lambda
    n <- length(y)
    meanx <- apply(X,2,mean)
    normx <- sqrt(apply((t(X)-meanx)^2,1,sum)/n)
    nz <- which(normx > .0001)
    XX <- scale(X[,nz],meanx[nz],normx[nz])
    p <- ncol(XX)
    if (family=="gaussian") yy <- y - mean(y)
    else yy <- y
    if (missing(lambda))
      {
        lambda <- setupLambda(XX,yy,family,alpha,lambda.min,nlambda)
        user.lambda <- FALSE
      }
    else
      {
        nlambda <- length(lambda)
        user.lambda <- TRUE
      }
    rm(XX)

    error <- array(NA,dim=c(nfolds,length(lambda)))
    
    if (family=="gaussian")
      {
        cv.ind <- ceiling((1:n)/n*nfolds)
      }
    else if (family=="binomial")
      {
        ind1 <- which(y==1)
        ind0 <- which(y==0)
        n1 <- length(ind1)
        n0 <- length(ind0)
        cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
        cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
        cv.ind <- numeric(n)
        cv.ind[y==1] <- cv.ind1
        cv.ind[y==0] <- cv.ind0
      }

    for (i in 1:nfolds)
      {
        if (trace) cat("Starting CV fold #",i,sep="","\n")
        X1 <- X[cv.ind!=i,]
        y1 <- y[cv.ind!=i]
        X2 <- X[cv.ind==i,]
        y2 <- y[cv.ind==i]

        fit.i <- ncvreg(X1,y1,family=family,alpha=alpha,lambda=lambda,warn=FALSE,...)
        yhat <- predict(fit.i,X2,type="response")
        error[i,1:ncol(yhat)] <- loss.ncvreg(y2,yhat,family)
      }

    ## Eliminate saturated lambda values, if any
    ind <- which(apply(is.finite(error),2,all))
    E <- error[,ind]
    lambda <- lambda[ind]
    
    val <- list(E=E,
                cve=apply(E,2,mean),
                lambda=lambda)
    val$min <- which.min(val$cve)
    class(val) <- "cv.ncvreg"
    return(val)
  }

