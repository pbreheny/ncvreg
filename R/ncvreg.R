ncvreg <- function(X, y, family=c("gaussian","binomial"), penalty=c("MCP", "SCAD", "lasso"), gamma=3, alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=.001, max.iter=1000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(X)), warn=TRUE, ...)
{
  ## Error checking
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  
  ## Deprication support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda
  
  ## Set up XX, yy, lambda
  XX <- standardize(X)
  center <- attr(XX, "center")
  scale <- attr(XX, "scale")
  nz <- which(scale > 1e-6)
  XX <- XX[ ,nz, drop=FALSE]
  yy <- if (family=="gaussian") y - mean(y) else y
  p <- ncol(XX)
  n <- length(yy)
  penalty.factor <- penalty.factor[nz]
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, family, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  if (family=="gaussian") {
    res <- .Call("cdfit_gaussian", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    b <- rbind(mean(y), matrix(res[[1]], p, nlambda))
    loss <- res[[2]]
    iter <- res[[3]]
  } else if (family=="binomial") {
    res <- .Call("cdfit_binomial", XX, as.numeric(yy), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    b <- rbind(res[[1]], matrix(res[[2]], p, nlambda))
    loss <- res[[3]]
    iter <- res[[4]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(b[p,])
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), family) else NULL
  
  ## Unstandardize
  b <- unstandardize(b, center[nz], scale[nz])
  beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
  beta[1,] <- b[1,]
  beta[nz+1,] <- b[-1,]
  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  varnames <- c("(Intercept)", varnames)
  dimnames(beta) <- list(varnames, round(lambda,digits=4))
  
  ## Output
  structure(list(beta = beta,
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
}
