ncvsurv <- function(X, y, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, 3),
                    alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=.001, max.iter=1000,
                    convex=TRUE, dfmax=p, penalty.factor=rep(1, ncol(X)), warn=TRUE, returnX=FALSE, ...) {

  # Coersion
  penalty <- match.arg(penalty)
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (class(y) != "matrix") {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must be a matrix or able to be coerced to a matrix")
    if (ncol(y)!=2) stop("y must have two columns for survival data: time-on-study and a censoring indicator")
  }
  if (storage.mode(y)=="integer") storage.mode(y) <- "double"
  if (storage.mode(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"

  ## Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")

  ## Set up XX, yy, lambda
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]
  ind <- order(y[,1])
  yy <- as.numeric(y[ind,1])
  Delta <- y[ind,2]
  XX <- XX[ind,,drop=FALSE]
  n <- length(yy)
  p <- ncol(XX)
  penalty.factor <- penalty.factor[nz]
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XX, yy, Delta, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  res <- .Call("cdfit_cox_dh", XX, yy, Delta, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor,
               alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
  b <- matrix(res[[1]], p, nlambda)
  loss <- -1*res[[2]]
  iter <- res[[3]]

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for some values of lambda")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), "cox", penalty.factor, Delta=Delta) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
  bb <- b/scale[nz]
  beta[nz,] <- bb
  offset <- -crossprod(center[nz], bb)

  ## Names
  varnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  dimnames(beta) <- list(varnames, round(lambda,digits=4))

  ## Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n),
                   class = c("ncvsurv", "ncvreg"))
  val$W <- exp(sweep(XX %*% b, 2, offset, "-"))
  val$time <- yy
  val$fail <- Delta
  if (returnX) {
    val$X <- XX
    val$center <- center
    val$scale <- scale
    val$y <- yy
  }
  val
}
