cv.ncvreg <- function(X, y, ..., cluster, nfolds=10, seed, fold, returnY=FALSE, trace=FALSE) {

  # Coercion
  if ('cv.ind' %in% names(list(...))) {
    seed <- list(...)$cv.ind
    warning("cv.ind has been deprecated and renamed fold; please use fold= in the future,\n
as support for cv.ind() is likely to be discontinued at some point.")
  }
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }

  fit <- ncvreg(X=X, y=y, ...)
  n <- length(y)
  E <- Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))
  if (fit$family=="binomial") {
    PE <- E
    if (!identical(sort(unique(y)), 0:1)) y <- as.numeric(y==max(y))
  }

  if (!missing(seed)) set.seed(seed)
  sde <- sqrt(.Machine$double.eps)
  if (missing(fold)) {
    if (fit$family=="binomial") {
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      fold1 <- ceiling(sample(1:n1)/(n1+sde)*nfolds)
      fold0 <- ceiling(sample(1:n0)/(n0+sde)*nfolds)
      fold <- numeric(n)
      fold[y==1] <- fold1
      fold[y==0] <- fold0
    } else {
      fold <- ceiling(sample(1:n)/(n+sqrt(.Machine$double.eps))*nfolds)
    }
  }

  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  if (!missing(cluster)) {
    if (!("cluster" %in% class(cluster))) stop("cluster is not of class 'cluster'; see ?makeCluster")
    parallel::clusterExport(cluster, c("fold","fit","X", "y", "cv.args"), envir=environment())
    parallel::clusterCall(cluster, function() require(ncvreg))
    fold.results <- parallel::parLapply(cl=cluster, X=1:nfolds, fun=cvf, XX=X, y=y, fold=fold, cv.args=cv.args)
  }

  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
    } else {
      if (trace) cat("Starting CV fold #",i,sep="","\n")
      res <- cvf(i, X, y, fold, cv.args)
    }
    E[fold==i, 1:res$nl] <- res$loss
    if (fit$family=="binomial") PE[fold==i, 1:res$nl] <- res$pe
    Y[fold==i, 1:res$nl] <- res$yhat
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[,ind]
  lambda <- fit$lambda[ind]

  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)

  # Bias correction
  e <- sapply(1:nfolds, function(i) apply(E[fold==i,,drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  val <- list(cve=cve, cvse=cvse, fold=fold, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min],
              null.dev=mean(loss.ncvreg(y, rep(mean(y), n), fit$family)), Bias=Bias)
  if (fit$family=="binomial") {
    pe <- apply(PE, 2, mean)
    val$pe <- pe[is.finite(pe)]
  }
  if (returnY) val$Y <- Y
  structure(val, class="cv.ncvreg")
}
cvf <- function(i, XX, y, fold, cv.args) {
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  fit.i <- do.call("ncvreg", cv.args)

  X2 <- XX[fold==i, , drop=FALSE]
  y2 <- y[fold==i]
  yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
  loss <- loss.ncvreg(y2, yhat, fit.i$family)
  pe <- if (fit.i$family=="binomial") {(yhat < 0.5) == y2} else NULL
  list(loss=loss, pe=pe, nl=length(fit.i$lambda), yhat=yhat)
}
