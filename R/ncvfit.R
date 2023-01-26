#' Direct interface for nonconvex penalized regression (non-pathwise)
#' 
#' This function is intended for users who know exactly what they're doing and want complete
#' control over the fitting process:
#' no standardization is applied, no intercept is included, no path is fit.
#' All of these things are best practices for data analysis, so if you are choosing not to do
#' them, you are on your own -- there is no guarantee that your results will be meaningful.
#' Some things in particular that you should pay attention to:
#'   * If your model has an intercept, it is up to you to (un)penalize it properly, typically
#'     by settings its corresponding element of `penalty.factor` to zero.
#'   * You should provide initial values for the coefficients; in nonconvex optimization,
#'     initial values are very important in determining which local solution an algorithm
#'     converges to.
#' 
#' At the moment, this function only works for least-squares loss functions.  Additional
#' functionality for other loss functions (logistic, Cox) is in development.
#' 
#' @param X                Design matrix; no intercept will be added, no standardization
#'                         will occur (n x p matrix)
#' @param y                Response vector (length n vector)
#' @param init             Initial values for beta.  Default: zero (length p vector)
#' @param r                Residuals corresponding to `init`; these will be calculated if not
#'                         supplied, but if they have already been calculated elsewhere, it is
#'                         more efficient to pass them as an argument. WARNING: If you supply
#'                         an incorrect value of `r`, the solution will be incorrect. (length
#'                         n vector)
#' @param xtx              X scales: the jth element should equal `crossprod(X[,j])/n`. These
#'                         will be calculated if not supplied, but if they have already been
#'                         calculated elsewhere, it is more efficient to pass them as an
#'                         argument.  In particular, if X is standardized, one should pass
#'                         `xtx = rep(1, p)`.  WARNING: If you supply an incorrect value of
#'                         `xtx`, the solution will be incorrect. (length p vector)
#' @param penalty          Penalty function to be applied, either "MCP" (default), "SCAD", or
#'                         "lasso")
#' @param gamma            Tuning parameter of the MCP/SCAD penalty, as in `ncvreg()`; default
#'                         is 3 for MCP and 3.7 for SCAD.
#' @param alpha            Tuning paramter controlling the ridge component of penalty, as in
#'                         `ncvreg()`; default is 1 (meaning no ridge penalty)
#' @param lambda           Regularization parameter value at which to estimate beta; must be
#'                         scalar -- for pathwise optimization, see `ncvreg()`
#' @param eps              Convergence threshhold. The algorithm iterates until the RMSD for
#'                         the change in linear predictors for each coefficient is less than
#'                         eps. Default is 1e-4.
#' @param max.iter         Maximum number of allowed iterations; if this number is reached,
#'                         algorithm will terminate prior to convergence.  Default: 1000.
#' @param penalty.factor   Multiplicative factor for the penalty applied to each coefficient,
#'                         as in `ncvreg()`.  In particular, note that if you include an
#'                         intercept, you probably want to set its entry to zero here.
#' @param warn             Return warning messages for failures to converge and model
#'                         saturation? Default is TRUE.
#' 
#' @return A list containing:
#' * `beta`: The estimated regression coefficients
#' * `iter`: The number of iterations required to solve for `beta
#' * `loss`: The loss (residual sum of squares) at convergence
#' * `resid`: The residuals at convergence
#' * `lambda`: See above
#' * `penalty`: See above
#' * `gamma`: See above
#' * `alpha`: See above
#' * `penalty.factor`: See above
#' * `n`: Sample size
#' 
#' @examples 
#' data(Prostate)
#' X <- cbind(1, Prostate$X)
#' y <- Prostate$y
#' fit <- ncvfit(X, y, lambda=0.1, penalty.factor=c(0, rep(1, ncol(X)-1)))
#' fit$beta
#' # Compare with:
#' coef(ncvreg(X, y), 0.1)
#' # The unstandardized version makes little sense here, as it fails to account
#' # for differences in the scales of the predictors.
#' @export ncvfit
ncvfit <- function(X, y, init=rep(0, ncol(X)), r, xtx, penalty=c("MCP", "SCAD", "lasso"),
                    gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda, eps=1e-5,
                    max.iter=1000, penalty.factor=rep(1, ncol(X)), warn=TRUE) {
  # Coersion
  # family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (typeof(X)=="integer") storage.mode(X) <- "double"
  if (typeof(X)=="character") stop("X must be a numeric matrix", call.=FALSE)
  if (!is.null(ncol(y)) && ncol(y) > 1) stop("y should be a vector of responses, not a matrix", call.=FALSE)
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }
  if (!is.double(penalty.factor)) penalty.factor <- as.double(penalty.factor)
  if (nrow(X) != length(y)) stop("y and X must have the same number of observations", call.=FALSE)
  if (!missing(r) && length(r) != nrow(X)) stop("r must be a length n vector", call.=FALSE)
  if (!missing(xtx) && length(xtx) != ncol(X)) stop("xtx must be a length p vector", call.=FALSE)
  
  # Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  # if (family=="binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
  # if (family=="binomial" & !identical(sort(unique(y)), 0:1)) y <- as.double(y==max(y))
  if (length(lambda) != 1) stop("lambda must be length 1")

  # Handle r, xtx
  
  # Fit
  if (missing(r)) r <- rep(NA_real_, nrow(X))
  if (!is.double(r)) r <- as.double(r)
  if (missing(xtx)) xtx <- rep(NA_real_, ncol(X))
  if (!is.double(xtx)) xtx <- as.double(xtx)
  res <- .Call("rawfit_gaussian", X, y, init, r, xtx, penalty, as.double(lambda), eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha)
  beta <- res[[1]]
  loss <- res[[2]]
  iter <- res[[3]]
  resid <- res[[4]]
  if (warn & (iter==max.iter)) warning("Maximum number of iterations reached")
  
  # Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  names(beta) <- varnames
  
  # Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        resid = resid,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = length(y)))
  #if (family=="poisson") val$y <- y
  #if (family=="binomial") val$Eta <- Eta
  val
}
