ncvreg_raw <- function(X, y, family=c("gaussian","binomial","poisson"), penalty=c("MCP", "SCAD", "lasso"),
                   gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda.min=ifelse(n > p,.001,.05), nlambda=100,
                   lambda, eps=1e-4, max.iter=10000, convex=TRUE, dfmax= p+1, penalty.factor=rep(1, ncol(X)),
                   warn=TRUE, returnX=FALSE, intercept=TRUE, ...) {
  # Coersion
  family <- match.arg(family)
  if (family != "gaussian") stop("ncvreg_raw currently only support gaussian family (least squares)")
  penalty <- match.arg(penalty)
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }
  if (storage.mode(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"

  # Error checking
  if (gamma <= 1 && penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 && penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (any(is.na(y)) || any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  #if (family=="binomial" && length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data")
  #if (family=="binomial" && !identical(sort(unique(y)), 0:1)) y <- as.numeric(y==max(y))
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations")

  ## Deprication support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

  ## Set up XX, yy, lambda
  col.se <- apply(X, 2, sd)
  ns <- which(col.se > 1e-6)
  XX <- X[, ns, drop=FALSE]
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)

  if (intercept == TRUE) {
    XX <- cbind(1, XX)
    penalty.factor <- c(0, penalty.factor)
    ncoef <- p+1
  } else {
    ncoef <- p
  }

  if (family=="gaussian" && intercept == TRUE) {
    yy <- y - mean(y)
  } else {
    yy <- y
  }
  n <- length(yy)
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, family, alpha, lambda.min, nlambda, penalty.factor, raw=TRUE)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  #if (family=="gaussian") {
    res <- .Call("cdfit_raw", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    beta <- matrix(res[[1]], ncoef, nlambda)
    loss <- res[[2]]
    iter <- res[[3]]    
  #}

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  beta <- beta[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  #if (family=="binomial") Eta <- Eta[,ind]
  if (warn && sum(iter) == max.iter) warning("Maximum number of iterations reached")

  ## Local convexity?
  convex.min <- if (convex) convexMin(beta, XX, penalty, gamma, lambda*(1-alpha), family, penalty.factor, a=a) else NULL

  # Names
  varnames <- colnames(X)
  if (intercept == TRUE) {
    if (is.null(varnames)) varnames <- paste0("V",seq(p-1))
    varnames <- c("(Intercept)", varnames)
  } else if (is.null(varnames)) {
    varnames <- paste0("V",seq(p))
  }
  dimnames(beta) <- list(varnames, lamNames(lambda))

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
  #if (family=="poisson") val$y <- y
  #if (family=="binomial") val$Eta <- Eta
  if (returnX) {
    val$X <- XX
    val$y <- yy
  }
  val
}
