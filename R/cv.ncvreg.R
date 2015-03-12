cv.ncvreg <- function(X, y, ..., nfolds=10, seed, cv.ind, trace=FALSE) {

  ## Error checking
  if (!missing(seed)) set.seed(seed)
  if (class(X) != "matrix") {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }  
  
  fit <- ncvreg(X=X, y=y, ...)
  n <- length(y)
  E <- matrix(NA, nrow=n, ncol=length(fit$lambda))
  if (fit$family=="binomial") {
    ##if (min(table(y)) < nfolds) stop("nfolds is larger than the smaller of 0/1 in the data; decrease nfolds")
    PE <- E
  }
  
  if (missing(cv.ind)) {
    if (fit$family=="binomial" & (min(table(y)) > nfolds)) {
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
      cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
      cv.ind <- numeric(n)
      cv.ind[y==1] <- cv.ind1
      cv.ind[y==0] <- cv.ind0
    } else {
      cv.ind <- ceiling(sample(1:n)/n*nfolds)
    }
  }
    
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")

    cv.args <- list(...)
    cv.args$X <- X[cv.ind!=i, , drop=FALSE]
    cv.args$y <- y[cv.ind!=i]
    cv.args$lambda <- fit$lambda
    cv.args$warn <- FALSE
    fit.i <- do.call("ncvreg", cv.args)

    X2 <- X[cv.ind==i, , drop=FALSE]
    y2 <- y[cv.ind==i]
    yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
    
    E[cv.ind==i, 1:length(fit.i$lambda)] <- loss.ncvreg(y2, yhat, fit$family)
    if (fit$family=="binomial") PE[cv.ind==i, 1:length(fit.i$lambda)] <- (yhat < 0.5) == y2
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[,ind]
  lambda <- fit$lambda

  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)
  
  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=mean(loss.ncvreg(y, rep(mean(y), n), fit$family)))
  if (fit$family=="binomial") {
    pe <- apply(PE, 2, mean)
    val$pe <- pe[is.finite(pe)]
  }
  structure(val, class="cv.ncvreg")
}
