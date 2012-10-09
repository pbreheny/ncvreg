ncvreg <- function(X, y, family=c("gaussian","binomial"), penalty=c("MCP","SCAD"), gamma=3, alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=.001, max.iter=1000, convex=TRUE, dfmax=p+1, warn=TRUE)
{
  ## Error checking
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  
  ## Set up XX, yy, lambda
  XX <- standardize(X)
  center <- attr(XX, "center")
  scale <- attr(XX, "scale")
  nz <- which(scale > 1e-6)
  XX <- XX[ ,nz, drop=FALSE]
  yy <- if (family=="gaussian") y - mean(y) else y
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, family, alpha, lambda.min, nlambda)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  ## Fit
  n <- length(yy)
  p <- ncol(XX)
  if (family=="gaussian") {
    fit <- .C("cdfit_gaussian", double(p*nlambda), double(nlambda), integer(nlambda), as.double(XX), as.double(yy), as.integer(n), as.integer(p), penalty, as.double(lambda), as.integer(nlambda), as.double(eps), as.integer(max.iter), as.double(gamma), as.double(alpha), as.integer(dfmax), as.integer(user.lambda))
    b <- rbind(mean(y), matrix(fit[[1]],nrow=p))
    loss <- fit[[2]]
    iter <- fit[[3]]
  } else if (family=="binomial") {
    fit <- .C("cdfit_binomial",double(nlambda),double(p*nlambda),double(nlambda),integer(nlambda),as.double(XX),as.double(yy),as.integer(n),as.integer(p),penalty,as.double(lambda),as.integer(nlambda),as.double(eps),as.integer(max.iter),as.double(gamma),as.double(alpha),as.integer(dfmax),as.integer(user.lambda),as.integer(warn))
    b <- rbind(fit[[1]],matrix(fit[[2]],nrow=p))
    loss <- fit[[3]]
    iter <- fit[[4]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(b[p,])
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")
  
  convex.min <- if (convex) convexMin(b,XX,penalty,gamma,lambda*(1-alpha),family) else NULL
  
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
                 convex.min = convex.min,
                 loss = loss,
                 n = n),
            class = "ncvreg")
}
