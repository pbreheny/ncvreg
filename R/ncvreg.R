ncvreg <- function(X, y, family=c("gaussian","binomial","poisson"), penalty=c("MCP", "SCAD", "lasso"),
                   gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100,
                   lambda, eps=1e-4, max.iter=10000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(X)),
                   warn=TRUE, returnX, ...) {
  # Coersion
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (typeof(X)=="integer") storage.mode(X) <- "double"
  if (typeof(X)=="character") stop("X must be a numeric matrix", call.=FALSE)
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }
  if (!is.double(penalty.factor)) penalty.factor <- as.double(penalty.factor)

  # Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  if (family=="binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
  if (family=="binomial" & !identical(sort(unique(y)), 0:1)) y <- as.double(y==max(y))
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)
  
  ## Deprication support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

  ## Set up XX, yy, lambda
  XX <- std(X)
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  if (family=="gaussian") {
    yy <- y - mean(y)
  } else {
    yy <- y
  }
  n <- length(yy)
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, family, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    if (nlambda == 1) warning(lambda.warning, call.=FALSE)
    if (!is.double(lambda)) lambda <- as.double(lambda)
    user.lambda <- TRUE
  }
  
  # Allow local_mfdr shortcut; probably not the ideal way to handle this
  if (sys.nframe() > 1) {
    cl <- sys.call(-1)[[1]]
    if (is.name(cl)) {
      if (cl == 'local_mfdr') return(list(X=XX, y=yy))
    }
  }

  ## Fit
  if (family=="gaussian") {
    res <- .Call("cdfit_gaussian", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    a <- rep(mean(y), nlambda)
    b <- matrix(res[[1]], p, nlambda)
    loss <- res[[2]]
    iter <- res[[3]]
  } else {
    if (family=="binomial") {
      res <- .Call("cdfit_binomial", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    } else if (family=="poisson") {
      res <- .Call("cdfit_poisson", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    }
    a <- res[[1]]
    b <- matrix(res[[2]], p, nlambda)
    loss <- res[[3]]
    Eta <- matrix(res[[4]], n, nlambda)
    iter <- res[[5]]
  }

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  a <- a[ind]
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (family=="binomial") Eta <- Eta[, ind]
  if (warn & sum(iter)==max.iter) warning("Maximum number of iterations reached")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), family, penalty.factor, a=a) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns+1,] <- bb
  beta[1,] <- a - crossprod(attr(XX, "center")[ns], bb)

  ## Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  varnames <- c("(Intercept)", varnames)
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
  if (family=="poisson") val$y <- y
  if (family=="binomial") val$Eta <- Eta
  if (missing(returnX)) {
    if (utils::object.size(XX) > 1e8) {
      warning("Due to the large size of X (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      returnX <- TRUE
    }
  }
  if (returnX) {
    val$X <- XX
    val$y <- yy
  }
  val
}

lambda.warning <- 'ncvreg() is intended for pathwise optimization, not for single values of lambda.
  1. You are strongly encouraged to fit a path and extract the solution at the lambda value of interest, rather than use ncvreg() in this way.
  2. In particular, if you are using the MCP or SCAD penalties, be aware that you greatly increase your risk of converging to an inferior local maximum if you do not fit an entire path.
  3. You may wish to look at the ncvfit() function, which is intended for non-path (i.e., single-lambda) optimization and allows the user to supply initial values.'
