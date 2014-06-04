ncvreg <- function(X, y, family=c("gaussian","binomial","poisson","cox"), penalty=c("MCP", "SCAD", "lasso"), 
                   gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, 
                   lambda, eps=.001, max.iter=1000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(X)), 
                   warn=TRUE, returnX=FALSE, ...) {
  ## Error checking
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")
  if (family=="binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data")
  
  ## Deprication support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda
  
  ## Set up XX, yy, lambda
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]
  if (family=="gaussian") {
    yy <- as.numeric(y - mean(y))
  } else if (family=="cox") {
    ind <- order(y[,1])
    yy <- as.numeric(y[ind,1])
    Delta <- y[ind,2]
    XX <- XX[ind,,drop=FALSE]
  } else {
    yy <- as.numeric(y)
  }
  n <- length(yy)
  p <- ncol(XX)
  penalty.factor <- penalty.factor[nz]
  if (missing(lambda)) {
    if (family!="cox") {
      lambda <- setupLambda(XX, yy, family, alpha, lambda.min, nlambda, penalty.factor)
    } else {
      lambda <- setupLambdaCox(XX, yy, Delta, alpha, lambda.min, nlambda, penalty.factor)
    }
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  if (family=="gaussian") {
    res <- .Call("cdfit_gaussian", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    a <- rep(mean(y),nlambda)
    b <- matrix(res[[1]], p, nlambda)
    loss <- res[[2]]
    iter <- res[[3]]
  } else if (family=="binomial") {
    res <- .Call("cdfit_binomial", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    a <- res[[1]]
    b <- matrix(res[[2]], p, nlambda)
    loss <- res[[3]]
    iter <- res[[4]]
  } else if (family=="poisson") {
    res <- .Call("cdfit_poisson", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    a <- res[[1]]
    b <- matrix(res[[2]], p, nlambda)
    loss <- res[[3]]
    iter <- res[[4]]
  } else if (family=="cox") {
    res <- .Call("cdfit_cox_dh", XX, yy, Delta, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, 
                 alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    b <- matrix(res[[1]], p, nlambda)
    loss <- res[[2]]
    iter <- res[[3]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  if (family !="cox") a <- a[ind]
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for some values of lambda")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), family, penalty.factor, a=a, Delta=Delta) else NULL
  
  ## Unstandardize
  if (family=="cox") {
    beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
    beta[nz,] <- b / scale[nz]
  } else {
    ##b <- unstandardize(a, b, center[nz], scale[nz])
    beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
    bb <- b/scale[nz]
    beta[nz+1,] <- bb
    beta[1,] <- a - crossprod(center[nz], bb)
  }
  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  varnames <- if (family=="cox") varnames else c("(Intercept)", varnames)
  dimnames(beta) <- list(varnames, round(lambda,digits=4))
  
  ## Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        family = family,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n),
                   class = "ncvreg")
  if (family=="poisson") val$y <- y
  if (returnX) {
    val$X <- XX
    val$center <- center
    val$scale <- scale
    val$y <- yy
  }
  val
}
