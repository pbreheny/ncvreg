#' Fit an MCP- or SCAD-penalized survival model
#'
#' Fit coefficients paths for MCP- or SCAD-penalized Cox regression models over
#' a grid of values for the regularization parameter lambda, with option for an
#' additional L2 penalty.
#'
#' The sequence of models indexed by the regularization parameter \code{lambda}
#' is fit using a coordinate descent algorithm.  In order to accomplish this,
#' the second derivative (Hessian) of the Cox partial log-likelihood is
#' diagonalized (see references for details).  The objective function is
#' defined to be \deqn{Q(\beta|X, y) = \frac{1}{n} L(\beta|X, y) + }{Q(\beta|X,
#' y) = (1/n)*L(\beta|X, y) + P(\beta, \lambda),}\deqn{
#' P_\lambda(\beta)}{Q(\beta|X, y) = (1/n)*L(\beta|X, y) + P(\beta, \lambda),}
#' where the loss function L is the deviance (-2 times the partial
#' log-likelihood) from the Cox regression mode. See
#' [here](https://pbreheny.github.io/ncvreg/articles/web/models.html) for more
#' details.
#'
#' Presently, ties are not handled by \code{ncvsurv} in a particularly
#' sophisticated manner.  This will be improved upon in a future release of
#' \code{ncvreg}.
#'
#' @param X The design matrix of predictor values.  \code{ncvsurv} standardizes
#'   the data prior to fitting.
#' @param y The time-to-event outcome, as a two-column matrix or
#'   \code{\link[survival]{Surv}} object.  The first column should be time on
#'   study (follow up time); the second column should be a binary variable with
#'   1 indicating that the event has occurred and 0 indicating (right)
#'   censoring.
#' @param penalty The penalty to be applied to the model.  Either "MCP" (the
#'   default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details).
#'   Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the
#'   relative contributions from the MCP/SCAD penalty and the ridge, or L2
#'   penalty.  \code{alpha=1} is equivalent to MCP/SCAD penalty, while
#'   \code{alpha=0} would be equivalent to ridge regression.  However,
#'   \code{alpha=0} is not supported; \code{alpha} may be arbitrarily small, but
#'   not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max.
#'   Default is .001 if the number of observations is larger than the number of
#'   covariates and .05 otherwise.
#' @param nlambda The number of lambda values.  Default is 100.
#' @param lambda A user-specified sequence of lambda values.  By default, a
#'   sequence of values of length \code{nlambda} is computed, equally spaced on
#'   the log scale.
#' @param eps Convergence threshhold.  The algorithm iterates until the RMSD for
#'   the change in linear predictors for any coefficient is less than
#'   \code{eps}.  Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path).
#'   Default is 1000.
#' @param convex Calculate index for which objective function ceases to be
#'   locally convex?  Default is TRUE.
#' @param dfmax Upper bound for the number of nonzero coefficients.  Default is
#'   no upper bound.  However, for large data sets, computational burden may be
#'   heavy for models with a large number of nonzero coefficients.
#' @param penalty.factor A multiplicative factor for the penalty applied to each
#'   coefficient.  If supplied, \code{penalty.factor} must be a numeric vector
#'   of length equal to the number of columns of \code{X}.  The purpose of
#'   \code{penalty.factor} is to apply differential penalization if some
#'   coefficients are thought to be more likely than others to be in the model.
#'   In particular, \code{penalty.factor} can be 0, in which case the
#'   coefficient is always in the model without any penalization/shrinkage.
#' @param warn Return warning messages for failures to converge and model
#'   saturation?  Default is TRUE.
#' @param returnX Return the standardized design matrix along with the fit?  By
#'   default, this option is turned on if X is under 100 MB, but turned off for
#'   larger matrices to preserve memory.  Note that certain methods, such as
#'   [summary.ncvreg()], require access to the design matrix and may not be able
#'   to run if `returnX=FALSE`.
#' @param ... Not used.
#'
#' @return An object with S3 class `ncvsurv` containing:
#' \describe{
#' \item{beta}{The fitted matrix of coefficients.  The number of rows is equal
#' to the number of coefficients, and the number of columns is equal to
#' \code{nlambda}.} \item{iter}{A vector of length \code{nlambda} containing
#' the number of iterations until convergence at each value of \code{lambda}.}
#' \item{lambda}{The sequence of regularization parameter values in the path.}
#' \item{penalty}{Same as above.} \item{model}{Same as above.}
#' \item{gamma}{Same as above.} \item{alpha}{Same as above.}
#' \item{convex.min}{The last index for which the objective function is locally
#' convex.  The smallest value of lambda for which the objective function is
#' convex is therefore \code{lambda[convex.min]}, with corresponding
#' coefficients \code{beta[,convex.min]}.} \item{loss}{The deviance of the
#' fitted model at each value of \code{lambda}.} \item{penalty.factor}{Same as
#' above.} \item{n}{The number of observations.}
#' }
#'
#'   For Cox models, the following objects are also returned (and are necessary
#'   to estimate baseline survival conditonal on the estimated regression
#'   coefficients), all of which are ordered by time on study.  I.e., the ith
#'   row of \code{W} does not correspond to the ith row of \code{X}):
#'
#' \describe{
#' \item{W}{Matrix of `exp(beta)` values for each subject over all `lambda` values.}
#' \item{time}{Times on study.}
#' \item{fail}{Failure event indicator.}
#' }
#'
#'   Additionally, if `returnX=TRUE`, the object will also contain
#'
#' \describe{
#' \item{X}{The standardized design matrix.}
#' }
#'
#' @seealso [plot.ncvreg()], [cv.ncvsurv()]
#'
#' @references
#' \itemize{
#' \item Breheny P and Huang J. (2011) Coordinate descent algorithms for nonconvex
#' penalized regression, with applications to biological feature selection.
#' *Annals of Applied Statistics*, **5**: 232-253. \doi{10.1214/10-AOAS388}
#'
#' \item Simon N, Friedman JH, Hastie T, and Tibshirani R. (2011)
#' Regularization Paths for Cox's Proportional Hazards Model via Coordinate
#' Descent. *Journal of Statistical Software*, **39**: 1-13.
#' \doi{10.18637/jss.v039.i05}
#' }
#' 
#' @examples
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#'
#' op <- par(mfrow=c(2,2))
#' fit <- ncvsurv(X, y)
#' plot(fit, main=expression(paste(gamma,"=",3)))
#' fit <- ncvsurv(X, y, gamma=10)
#' plot(fit, main=expression(paste(gamma,"=",10)))
#' fit <- ncvsurv(X, y, gamma=1.5)
#' plot(fit, main=expression(paste(gamma,"=",1.5)))
#' fit <- ncvsurv(X, y, penalty="SCAD")
#' plot(fit, main=expression(paste("SCAD, ",gamma,"=",3)))
#' par(op)
#'
#' fit <- ncvsurv(X,y)
#' ll <- log(fit$lambda)
#' op <- par(mfrow=c(2,1))
#' plot(ll, BIC(fit), type="l", xlim=rev(range(ll)))
#' lam <- fit$lambda[which.min(BIC(fit))]
#' b <- coef(fit, lambda=lam)
#' b[b!=0]
#' plot(fit)
#' abline(v=lam)
#' par(op)
#'
#' S <- predict(fit, X, type='survival', lambda=lam)
#' plot(S, xlim=c(0,200))
#' @export

ncvsurv <- function(X, y, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, 3),
                    alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=1e-4, max.iter=10000,
                    convex=TRUE, dfmax=p, penalty.factor=rep(1, ncol(X)), warn=TRUE, returnX, ...) {

  # Coersion
  penalty <- match.arg(penalty)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (!inherits(y, "matrix")) {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be a matrix or able to be coerced to a matrix", call.=FALSE)
    if (ncol(y) != 2) stop("y must have two columns for survival data: time-on-study and a censoring indicator", call.=FALSE)
  }
  if (typeof(y) == "integer") storage.mode(y) <- "double"
  if (typeof(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"

  ## Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)

  ## Set up XX, yy, lambda
  tOrder <- order(y[, 1])
  yy <- as.double(y[tOrder, 1])
  Delta <- y[tOrder, 2]
  n <- length(yy)
  XX <- std(X[tOrder, , drop=FALSE])
  if (sys.nframe() > 1 && sys.call(-1)[[1]]=="local_mfdr") return(list(X=XX, time=yy, fail=Delta))
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XX, yy, Delta, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    if (!is.double(lambda)) lambda <- as.double(lambda)
    if (nlambda == 1) {
      warning(gsub('ncvreg', 'ncvsurv', lambda.warning), call.=FALSE)
    } else if (any(diff(lambda) > 0)) {
      lambda <- sort(lambda, decreasing=TRUE)
    }
    user.lambda <- TRUE
  }

  ## Fit
  res <- .Call("cdfit_cox_dh", XX, Delta, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor,
               alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
  b <- matrix(res[[1]], p, nlambda)
  loss <- -1*res[[2]]
  iter <- res[[3]]
  Eta <- matrix(res[[4]], n, nlambda)

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  Eta <- Eta[, ind, drop=FALSE]
  if (warn & sum(iter)==max.iter) warning("Algorithm failed to converge for some values of lambda")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), "cox", penalty.factor, Delta=Delta) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns,] <- bb
  offset <- -crossprod(attr(XX, "center")[ns], bb)

  ## Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  dimnames(beta) <- list(varnames, lam_names(lambda))
  obsnames <- if (is.null(rownames(X))) 1:nrow(X) else rownames(X)
  dimnames(Eta) <- list(obsnames, lam_names(lambda))
  
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
                        n = n,
                        time = yy,
                        fail = Delta,
                        order = tOrder,
                        linear.predictors = sweep(Eta, 2, colMeans(Eta), '-')),
                   class = c("ncvsurv", "ncvreg"))
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
  }
  val
}
