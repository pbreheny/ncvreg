cv.ncvreg <- function(X, y, ..., nfolds=10, seed, trace=FALSE) {
  if (!missing(seed)) set.seed(seed)
  fit <- ncvreg(X=X, y=y, ...)
  E <- matrix(NA, nrow=length(y), ncol=length(fit$lambda))
  if (fit$family=="binomial") {
    if (min(table(y)) < nfolds) stop("nfolds is larger than the smaller of 0/1 in the data; decrease nfolds")
    PE <- E
  }

  n <- length(y)
  if (fit$family=="gaussian" | fit$family=="poisson") {
    cv.ind <- ceiling(sample(1:n)/n*nfolds)
  } else if (fit$family=="binomial") {
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
    yhat <- predict(fit.i, X2, type="response")
    E[cv.ind==i, 1:ncol(yhat)] <- loss.ncvreg(y2, yhat, fit$family)
    if (fit$family=="binomial") PE[cv.ind==i, 1:ncol(yhat)] <- (yhat < 0.5) == y2
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
