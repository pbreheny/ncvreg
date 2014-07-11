cv.ncvsurv <- function(X, y, ..., nfolds=10, seed, trace=FALSE) {
  if (!missing(seed)) set.seed(seed)
  fit <- ncvsurv(X=X, y=y, ...)
  n <- nrow(X)
  E <- matrix(NA, nrow=n, ncol=length(fit$lambda))
  
  cv.ind <- ceiling(sample(1:n)/n*nfolds)

  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")

    cv.args <- list(...)
    cv.args$X <- X[cv.ind!=i, , drop=FALSE]
    cv.args$y <- y[cv.ind!=i,]
    cv.args$lambda <- fit$lambda
    cv.args$warn <- FALSE
    fit.i <- do.call("ncvsurv", cv.args)

    X2 <- X[cv.ind==i, , drop=FALSE]
    y2 <- y[cv.ind==i,]
    yhat <- predict(fit.i, X2, type="response")
    if (fit$model=="cox") {
      E[cv.ind==i, 1:ncol(yhat)] <- coxCVL(y, cv.ind, yhat)
    }
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[,ind]
  lambda <- fit$lambda

  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)
  
  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=NA)
  structure(val, class="cv.ncvsurv")
}
