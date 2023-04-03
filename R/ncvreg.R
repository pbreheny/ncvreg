#' Fit an MCP- or SCAD-penalized regression path
#' 
#' Fit coefficients paths for MCP- or SCAD-penalized regression models over a
#' grid of values for the regularization parameter lambda.  Fits linear and
#' logistic regression models, with option for an additional L2 penalty.
#' 
#' The sequence of models indexed by the regularization parameter \code{lambda}
#' is fit using a coordinate descent algorithm.  For logistic regression
#' models, some care is taken to avoid model saturation; the algorithm may exit
#' early in this setting.  The objective function is defined to be
#' \deqn{Q(\beta|X, y) = \frac{1}{n} L(\beta|X, y) + }{Q(\beta|X, y) =
#' (1/n)*L(\beta|X, y) + P(\beta, \lambda),}\deqn{ P_\lambda(\beta)}{Q(\beta|X,
#' y) = (1/n)*L(\beta|X, y) + P(\beta, \lambda),} where the loss function L is
#' the deviance (-2 times the log likelihood) for the specified outcome
#' distribution (gaussian/binomial/poisson). See
#' [here](https://pbreheny.github.io/ncvreg/articles/web/models.html) for more
#' details.
#' 
#' This algorithm is stable, very efficient, and generally converges quite
#' rapidly to the solution.  For GLMs,
#' [adaptive rescaling](https://myweb.uiowa.edu/pbreheny/pdf/Breheny2011.pdf)
#' is used.
#' 
#' @param X The design matrix, without an intercept.  \code{ncvreg}
#' standardizes the data and includes an intercept by default.
#' @param y The response vector.
#' @param family Either "gaussian", "binomial", or "poisson", depending on the
#' response.
#' @param penalty The penalty to be applied to the model.  Either "MCP" (the
#' default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details).
#' Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the
#' relative contributions from the MCP/SCAD penalty and the ridge, or L2
#' penalty.  \code{alpha=1} is equivalent to MCP/SCAD penalty, while
#' \code{alpha=0} would be equivalent to ridge regression.  However,
#' \code{alpha=0} is not supported; \code{alpha} may be arbitrarily small, but
#' not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of
#' lambda.max.  Default is .001 if the number of observations is larger than
#' the number of covariates and .05 otherwise.
#' @param nlambda The number of lambda values.  Default is 100.
#' @param lambda A user-specified sequence of lambda values.  By default, a
#' sequence of values of length \code{nlambda} is computed, equally spaced on
#' the log scale.
#' @param eps Convergence threshhold.  The algorithm iterates until the RMSD
#' for the change in linear predictors for each coefficient is less than
#' \code{eps}.  Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path).
#' Default is 10000.
#' @param convex Calculate index for which objective function ceases to be
#' locally convex?  Default is TRUE.
#' @param dfmax Upper bound for the number of nonzero coefficients.  Default is
#' no upper bound.  However, for large data sets, computational burden may be
#' heavy for models with a large number of nonzero coefficients.
#' @param penalty.factor A multiplicative factor for the penalty applied to
#' each coefficient.  If supplied, \code{penalty.factor} must be a numeric
#' vector of length equal to the number of columns of \code{X}.  The purpose of
#' \code{penalty.factor} is to apply differential penalization if some
#' coefficients are thought to be more likely than others to be in the model.
#' In particular, \code{penalty.factor} can be 0, in which case the coefficient
#' is always in the model without shrinkage.
#' @param warn Return warning messages for failures to converge and model
#' saturation?  Default is TRUE.
#' @param returnX Return the standardized design matrix along with the fit?  By
#' default, this option is turned on if X is under 100 MB, but turned off for
#' larger matrices to preserve memory.  Note that certain methods, such as
#' \code{\link{summary.ncvreg}} require access to the design matrix and may not
#' be able to run if \code{returnX=FALSE}.
#' @param ... Not used.
#' @return An object with S3 class \code{"ncvreg"} containing: \describe{
#' \item{beta}{The fitted matrix of coefficients.  The number of rows is equal
#' to the number of coefficients, and the number of columns is equal to
#' \code{nlambda}.} \item{iter}{A vector of length \code{nlambda} containing
#' the number of iterations until convergence at each value of \code{lambda}.}
#' \item{lambda}{The sequence of regularization parameter values in the path.}
#' \item{penalty}{Same as above.} \item{family}{Same as above.}
#' \item{gamma}{Same as above.} \item{alpha}{Same as above.}
#' \item{convex.min}{The last index for which the objective function is locally
#' convex.  The smallest value of lambda for which the objective function is
#' convex is therefore \code{lambda[convex.min]}, with corresponding
#' coefficients \code{beta[,convex.min]}.} \item{loss}{A vector containing the
#' deviance (i.e., the loss) at each value of \code{lambda}.  Note that for
#' \code{gaussian} models, the loss is simply the residual sum of squares.}
#' \item{penalty.factor}{Same as above.} \item{n}{Sample size.} }
#' 
#' Additionally, if \code{returnX=TRUE}, the object will also contain
#' 
#' \describe{ \item{X}{The standardized design matrix.} \item{y}{The response,
#' centered if \code{family='gaussian'}.} }
#' 
#' @seealso [plot.ncvreg()], [cv.ncvreg()]
#' 
#' @references Breheny P and Huang J. (2011) Coordinate descentalgorithms for
#' nonconvex penalized regression, with applications to biological feature
#' selection.  *Annals of Applied Statistics*, **5**: 232-253.
#' \doi{10.1214/10-AOAS388}
#'
#' @examples 
#' # Linear regression --------------------------------------------------
#' data(Prostate)
#' X <- Prostate$X
#' y <- Prostate$y
#' 
#' op <- par(mfrow=c(2,2))
#' fit <- ncvreg(X, y)
#' plot(fit, main=expression(paste(gamma,"=",3)))
#' fit <- ncvreg(X, y, gamma=10)
#' plot(fit, main=expression(paste(gamma,"=",10)))
#' fit <- ncvreg(X, y, gamma=1.5)
#' plot(fit, main=expression(paste(gamma,"=",1.5)))
#' fit <- ncvreg(X, y, penalty="SCAD")
#' plot(fit, main=expression(paste("SCAD, ",gamma,"=",3)))
#' par(op)
#' 
#' op <- par(mfrow=c(2,2))
#' fit <- ncvreg(X, y)
#' plot(fit, main=expression(paste(alpha,"=",1)))
#' fit <- ncvreg(X, y, alpha=0.9)
#' plot(fit, main=expression(paste(alpha,"=",0.9)))
#' fit <- ncvreg(X, y, alpha=0.5)
#' plot(fit, main=expression(paste(alpha,"=",0.5)))
#' fit <- ncvreg(X, y, alpha=0.1)
#' plot(fit, main=expression(paste(alpha,"=",0.1)))
#' par(op)
#' 
#' op <- par(mfrow=c(2,2))
#' fit <- ncvreg(X, y)
#' plot(mfdr(fit))             # Independence approximation
#' plot(mfdr(fit), type="EF")  # Independence approximation
#' perm.fit <- perm.ncvreg(X, y)
#' plot(perm.fit)
#' plot(perm.fit, type="EF")
#' par(op)
#' 
#' # Logistic regression ------------------------------------------------
#' data(Heart)
#' X <- Heart$X
#' y <- Heart$y
#' 
#' op <- par(mfrow=c(2,2))
#' fit <- ncvreg(X, y, family="binomial")
#' plot(fit, main=expression(paste(gamma,"=",3)))
#' fit <- ncvreg(X, y, family="binomial", gamma=10)
#' plot(fit, main=expression(paste(gamma,"=",10)))
#' fit <- ncvreg(X, y, family="binomial", gamma=1.5)
#' plot(fit, main=expression(paste(gamma,"=",1.5)))
#' fit <- ncvreg(X, y, family="binomial", penalty="SCAD")
#' plot(fit, main=expression(paste("SCAD, ",gamma,"=",3)))
#' par(op)
#' 
#' op <- par(mfrow=c(2,2))
#' fit <- ncvreg(X, y, family="binomial")
#' plot(fit, main=expression(paste(alpha,"=",1)))
#' fit <- ncvreg(X, y, family="binomial", alpha=0.9)
#' plot(fit, main=expression(paste(alpha,"=",0.9)))
#' fit <- ncvreg(X, y, family="binomial", alpha=0.5)
#' plot(fit, main=expression(paste(alpha,"=",0.5)))
#' fit <- ncvreg(X, y, family="binomial", alpha=0.1)
#' plot(fit, main=expression(paste(alpha,"=",0.1)))
#' par(op)
#' 
#' @export ncvreg

ncvreg <- function(X, y, family=c("gaussian","binomial","poisson"), penalty=c("MCP", "SCAD", "lasso"),
                   gamma=switch(penalty, SCAD=3.7, 3), alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100,
                   lambda, eps=1e-4, max.iter=10000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(X)),
                   warn=TRUE, returnX, ...) {
  # Coersion
  family <- match.arg(family)
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
    Eta <- matrix(res[[3]], nrow=n) + mean(y)
    iter <- res[[4]]
  } else {
    res <- .Call("cdfit_glm", XX, yy, family, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
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
  Eta <- Eta[, ind, drop=FALSE]
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
  dimnames(beta) <- list(varnames, lam_names(lambda))
  obsnames <- if (is.null(rownames(X))) 1:nrow(X) else rownames(X)
  dimnames(Eta) <- list(obsnames, lam_names(lambda))

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
                        linear.predictors = Eta,
                        penalty.factor = penalty.factor,
                        n = n,
                        y = y),
                   class = "ncvreg")
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

lambda.warning <- 'ncvreg() is intended for pathwise optimization, not for single values of lambda.
  1. You are strongly encouraged to fit a path and extract the solution at the lambda value of interest, rather than use ncvreg() in this way.
  2. In particular, if you are using the MCP or SCAD penalties, be aware that you greatly increase your risk of converging to an inferior local maximum if you do not fit an entire path.
  3. You may wish to look at the ncvfit() function, which is intended for non-path (i.e., single-lambda) optimization and allows the user to supply initial values.'
