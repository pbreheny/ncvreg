ncvreg_fit <- function(X, y, family=c("gaussian","binomial","poisson"), penalty=c("MCP", "SCAD", "lasso"), gamma=3, alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=.001, max.iter=1000, dfmax=p+1, penalty.factor=rep(1, ncol(X)), warn=TRUE) {
  
  ## Error checking
  if (class(X) != "matrix") {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  n <- length(y)
  p <- ncol(X)
  
  ## Set up lambda
  if (missing(lambda)) {
    lambda <- setupLambda(X, y, family, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    user.lambda <- TRUE
  }
  nlam <- length(lambda)
  
  if (family=="gaussian") {
    res <- .Call("cdfit_raw", X, y, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    b <- matrix(res[[1]], p, nlam)
    loss <- res[[2]]
    iter <- res[[3]]
  } else if (family=="binomial") {
    res <- .Call("cdfit_binomial", X, as.numeric(y), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    b0 <- numeric(nlam)
    b <- rbind(res[[1]], matrix(res[[2]], p, nlam))
    loss <- res[[3]]
    iter <- res[[4]]
  } else if (family=="poisson") {
    res <- .Call("cdfit_poisson", X, as.numeric(y), penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
    b0 <- numeric(nlam)
    b <- rbind(res[[1]], matrix(res[[2]], p, nlam))
    loss <- res[[3]]
    iter <- res[[4]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for some values of lambda")

  ## Output
  val <- structure(list(beta = b,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        family = family,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = NULL,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n),
                   class = "ncvreg")
  if (family=="poisson") val$y <- y
  val
}
