cv.ncvreg <- function(X, y, ..., cluster, nfolds=10, seed, cv.ind, returnY=FALSE, trace=FALSE) {

  # Coersion
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
  if (missing(cv.ind)) {
    if (fit$family=="binomial" & (min(table(y)) > nfolds)) {
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      cv.ind1 <- ceiling(sample(1:n1)/(n1+sde)*nfolds)
      cv.ind0 <- ceiling(sample(1:n0)/(n0+sde)*nfolds)
      cv.ind <- numeric(n)
      cv.ind[y==1] <- cv.ind1
      cv.ind[y==0] <- cv.ind0
    } else {
      cv.ind <- ceiling(sample(1:n)/(n+sqrt(.Machine$double.eps))*nfolds)
    }
  }

  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  if (!missing(cluster)) {
    if (!("cluster" %in% class(cluster))) stop("cluster is not of class 'cluster'; see ?makeCluster")
    parallel::clusterExport(cluster, c("cv.ind","fit","X", "y", "cv.args"), envir=environment())
    parallel::clusterCall(cluster, function() require(ncvreg))
    fold.results <- parallel::parLapply(cl=cluster, X=1:nfolds, fun=cvf, XX=X, y=y, cv.ind=cv.ind, cv.args=cv.args)
  }

  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
    } else {
      if (trace) cat("Starting CV fold #",i,sep="","\n")
      res <- cvf(i, X, y, cv.ind, cv.args)
    }
    E[cv.ind==i, 1:res$nl] <- res$loss
    if (fit$family=="binomial") PE[cv.ind==i, 1:res$nl] <- res$pe
    Y[cv.ind==i, 1:res$nl] <- res$yhat
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
  e <- sapply(1:nfolds, function(i) apply(E[cv.ind==i,,drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min],
              null.dev=mean(loss.ncvreg(y, rep(mean(y), n), fit$family)), Bias=Bias)
  if (fit$family=="binomial") {
    pe <- apply(PE, 2, mean)
    val$pe <- pe[is.finite(pe)]
  }
  if (returnY) val$Y <- Y
  structure(val, class="cv.ncvreg")
}
cvf <- function(i, XX, y, cv.ind, cv.args) {
  cv.args$X <- XX[cv.ind!=i, , drop=FALSE]
  cv.args$y <- y[cv.ind!=i]
  fit.i <- do.call("ncvreg", cv.args)

  X2 <- XX[cv.ind==i, , drop=FALSE]
  y2 <- y[cv.ind==i]
  yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
  loss <- loss.ncvreg(y2, yhat, fit.i$family)
  pe <- if (fit.i$family=="binomial") {(yhat < 0.5) == y2} else NULL
  list(loss=loss, pe=pe, nl=length(fit.i$lambda), yhat=yhat)
}
