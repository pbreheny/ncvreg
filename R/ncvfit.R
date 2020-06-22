#' Direct interface for nonconvex penalized regression (non-pathwise)
#' 
#' This function is intended for users who know exactly what they're doing and want complete control over the fitting process:
#' no standardization is applied, no intercept is included, no path is fit.
#' All of these things are best practices for data analysis, so if you are choosing not to do them, you are on your own -- I cannot guarantee that your results will be meaningful.
#' Some things in particular that you should pay attention to:
#'   * If your model has an intercept, it is up to you to (un)penalize it properly, typically by settings its corresponding element of `penalty.factor` to zero.
#'   * You should probably provide an initial value; in nonconvex optimization, initial values are very important in determining which local solution an algorithm converges to.
#' 
#' To-do: Add family

ncvfit <- function(X, y, init=rep(0, ncol(X)), penalty=c("MCP", "SCAD", "lasso"),
                    gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda, eps=1e-4,
                    max.iter=1000, penalty.factor=rep(1, ncol(X)), warn=TRUE) {
  # Coersion
  # family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (typeof(X)=="integer") storage.mode(X) <- "double"
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  if (!is.double(penalty.factor)) penalty.factor <- as.double(penalty.factor)
  
  # Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X", call.=FALSE)
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  # if (family=="binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
  # if (family=="binomial" & !identical(sort(unique(y)), 0:1)) y <- as.double(y==max(y))
  if (length(lambda) != 1) stop("lambda must be length 1")
  
  # Fit
  n <- length(y)
  res <- .Call("rawfit_gaussian", X, y, init, penalty, as.double(lambda), eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha)
  beta <- res[[1]]
  loss <- res[[2]]
  iter <- res[[3]]
  if (warn & (iter==max.iter)) warning("Maximum number of iterations reached")
  
  # Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  names(beta) <- varnames
  
  # Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n))
  #if (family=="poisson") val$y <- y
  #if (family=="binomial") val$Eta <- Eta
  val
}

