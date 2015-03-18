cv.ncvsurv <- function(X, y, ..., nfolds=10, seed, trace=FALSE, events.only=TRUE) {

  ## Error checking
  if (class(X) != "matrix") {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (class(y) != "matrix") {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must be a matrix or able to be coerced to a matrix")
    if (ncol(y)!=2) stop("y must have two columns for survival data: time-on-study and a censoring indicator")
  }

  fit <- ncvsurv(X=X, y=y, ...)
  n <- nrow(X)
  E <- matrix(NA, nrow=n, ncol=length(fit$lambda))
  
  if (!missing(seed)) set.seed(seed)
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
    if (fit$model=="cox") {
      eta <- predict(fit.i, X)
      ll <- loss.ncvsurv(y, eta)
      for (ii in which(cv.ind==i)) {
        eta.ii <- predict(fit.i, X[-ii,])
        E[ii, 1:ncol(eta)] <- -2*(ll-loss.ncvsurv(y[-ii,], eta.ii))
      }
    }
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[,ind]
  lambda <- fit$lambda

  ## Return
  if (events.only) E <- E[y[,2]==1,]
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(nrow(E))
  min <- which.min(cve)
  
  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=cve[1])
  structure(val, class=c("cv.ncvsurv", "cv.ncvreg"))
}
