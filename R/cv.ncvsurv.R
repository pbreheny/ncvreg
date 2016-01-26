cv.ncvsurv <- function(X, y, ..., cluster, nfolds=10, seed, returnY=FALSE, trace=FALSE) {

  # Complete data fit
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$returnX <- TRUE
  fit <- do.call("ncvsurv", fit.args)

  # Get standardized X, y
  X <- fit$X
  y <- cbind(fit$time, fit$fail)
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) fit$X <- NULL

  # Set up folds
  n <- nrow(X)
  if (!missing(seed)) set.seed(seed)
  cv.ind <- ceiling(sample(1:n)/n*nfolds)
  Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))

  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  cv.args$penalty.factor <- fit$penalty.factor
  if (!missing(cluster)) {
    if (!("cluster" %in% class(cluster))) stop("cluster is not of class 'cluster'; see ?makeCluster")
    parallel::clusterExport(cluster, c("cv.ind","fit","X", "y", "cv.args"), envir=environment())
    parallel::clusterCall(cluster, function() require(ncvreg))
    fold.results <- parallel::parLapply(cl=cluster, X=1:nfolds, fun=cvf.surv, XX=X, y=y, cv.ind=cv.ind, cv.args=cv.args)
  }

  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
    } else {
      if (trace) cat("Starting CV fold #",i,sep="","\n")
      res <- cvf.surv(i, X, y, cv.ind, cv.args)
    }
    Y[cv.ind==i, 1:res$nl] <- res$yhat
    #E[cv.ind==i, 1:res$nl] <- res$loss
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(Y), 2, all))
  #E <- E[,ind]
  Y <- Y[,ind]
  lambda <- fit$lambda[ind]

  ## Return
  cve <- as.numeric(loss.ncvsurv(y, Y))
  #if (events.only) E <- E[y[,2]==1,]
  #cve <- apply(E, 2, mean)
  #cvse <- apply(E, 2, sd) / sqrt(nrow(E))
  min <- which.min(cve)

  val <- list(cve=cve, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=cve[1])
  if (returnY) val$Y <- Y
  structure(val, class=c("cv.ncvsurv", "cv.ncvreg"))
}
cvf.surv <- function(i, XX, y, cv.ind, cv.args) {
  cv.args$X <- XX[cv.ind!=i, , drop=FALSE]
  cv.args$y <- y[cv.ind!=i,]
  fit.i <- do.call("ncvsurv", cv.args)

  X2 <- XX[cv.ind==i, , drop=FALSE]
  y2 <- y[cv.ind==i,]
  nl <- length(fit.i$lambda)
  yhat <- predict(fit.i, X2)

#   eta <- predict(fit.i, XX)
#   ll <- loss.ncvsurv(y, eta)
#   eta.i <- predict(fit.i, XX[cv.ind!=i, , drop=FALSE])
#   ll.i <- loss.ncvsurv(y[cv.ind!=i,], eta.i)
#   loss <- matrix(ll.i, sum(cv.ind==i), nl, byrow=TRUE)

#   ind <- which(cv.ind==i)
#   n <- length(ind)
#   loss <- matrix(NA, n, nl)
#   for (j in 1:n) {
#     k <- ind[j]
#     eta.j <- predict(fit.i, XX[-k,])
#     loss[j,] <- 2*(ll-loss.ncvsurv(y[-k,], eta.j))
#   }
  list(nl=length(fit.i$lambda), yhat=yhat)#, loss=loss)
}
