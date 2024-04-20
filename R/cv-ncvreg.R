#' Cross-validation for ncvreg/ncvsurv
#' 
#' Performs k-fold cross validation for MCP- or SCAD-penalized regression
#' models over a grid of values for the regularization parameter lambda.
#' 
#' The function calls `ncvreg`/`ncvsurv` `nfolds` times, each
#' time leaving out 1/`nfolds` of the data.  The cross-validation error is
#' based on the deviance; [see here for more details](https://pbreheny.github.io/ncvreg/articles/web/models.html).
#' 
#' For `family="binomial"` models, the cross-validation fold assignments are
#' balanced across the 0/1 outcomes, so that each fold has the same proportion
#' of 0/1 outcomes (or as close to the same proportion as it is possible to
#' achieve if cases do not divide evenly).
#' 
#' For Cox models, `cv.ncvsurv()` uses the approach of calculating the full
#' Cox partial likelihood using the cross-validated set of linear predictors.
#' Other approaches to cross-validation for the Cox regression model have been
#' proposed in the literature; the strengths and weaknesses of the various
#' methods for penalized regression in the Cox model are the subject of current
#' research.  A simple approximation to the standard error is provided,
#' although an option to bootstrap the standard error (`se='bootstrap'`) is also
#' available.
#' 
#' @aliases cv.ncvreg cv.ncvsurv
#' 
#' @param X         The design matrix, without an intercept, as in [ncvreg()] or [ncvsurv()].
#' @param y         The response, as in [ncvreg()] or [ncvsurv()].
#' @param ...       Additional arguments to [ncvreg()] or [ncvsurv()].
#' @param cluster   `cv.ncvreg()` and `cv.ncvsurv()` can be run in parallel
#' across a cluster using the **parallel** package. The cluster must be set
#' up in advance using the [parallel::makeCluster()] function from that package.
#' The cluster must then be passed to `cv.ncvreg()` or `cv.ncvsurv()` (see example).
#' @param nfolds    The number of cross-validation folds.  Default is 10.
#' @param fold      Which fold each observation belongs to. By default the observations are randomly assigned.
#' @param seed      You may set the seed of the random number generator in order to obtain reproducible results.
#' @param returnY   Should `cv.ncvreg()`/`cv.ncvsurv()` return the linear predictors
#'   from the cross-validation folds?  Default is `FALSE`; if `TRUE`, this will
#'   return a matrix in which the element for row i, column j is the fitted
#'   value for observation i from the fold in which observation i was excluded
#'   from the fit, at the jth value of lambda. NOTE: For `cv.ncvsurv()`, the
#'   rows of `Y` are ordered by time on study, and therefore will not correspond
#'   to the original order of observations pased to `cv.ncvsurv()`.
#' @param trace     If set to `TRUE`, inform the user of progress by announcing
#'   the beginning of each CV fold. Default is `FALSE`.
#' @param se        For `cv.ncvsurv()`, the method by which the cross-valiation
#'   standard error (CVSE) is calculated. The 'quick' approach is based on a
#'   rough approximation, but can be calculated more or less instantly.  The
#'   'bootstrap' approach is more accurate, but requires additional computing time.
#'   
#' @returns An object with S3 class `cv.ncvreg` or `cv.ncvsurv` containing:
#' \describe{
#'   \item{cve}{The error for each value of `lambda`, averaged across the cross-
#'     validation folds.}
#'   \item{cvse}{The estimated standard error associated with each value of for `cve`.}
#'   \item{fold}{The fold assignments for cross-validation for each observation;
#'     note that for `cv.ncvsurv()`, these are in terms of the ordered observations,
#'     not the original observations.}
#'   \item{lambda}{The sequence of regularization parameter values along which
#'     the cross-validation error was calculated.}
#'   \item{fit}{The fitted [ncvreg()] or [ncvsurv()] object for the whole data.}
#'   \item{min}{The index of `lambda` corresponding to `lambda.min`.}
#'   \item{lambda.min}{The value of `lambda` with the minimum cross-validation error.}
#'   \item{null.dev}{The deviance for the intercept-only model. If you have supplied
#'     your own `lambda` sequence, this quantity may not be meaningful.}
#'   \item{Bias}{The estimated bias of the minimum cross-validation error, as in
#'     Tibshirani and Tibshirani (2009) \doi{10.1214/08-AOAS224}}
#'   \item{pe}{If `family="binomial"`, the cross-validation prediction error for
#'     each value of `lambda`.}
#'   \item{Y}{If `returnY=TRUE`, the matrix of cross-validated fitted values (see above).}
#' }
#' 
#' @author Patrick Breheny; Grant Brown helped with the parallelization support
#' 
#' @seealso [ncvreg()], [plot.cv.ncvreg()], [summary.cv.ncvreg()]
#' 
#' @references
#' Breheny P and Huang J. (2011) Coordinate descent algorithms for nonconvex
#' penalized regression, with applications to biological feature selection.
#' *Annals of Applied Statistics*, **5**: 232-253. \doi{10.1214/10-AOAS388}
#' 
#' @examples
#' data(Prostate)
#' 
#' cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
#' plot(cvfit)
#' summary(cvfit)
#' 
#' fit <- cvfit$fit
#' plot(fit)
#' beta <- fit$beta[,cvfit$min]
#' 
#' ## requires loading the parallel package
#' \dontrun{
#' library(parallel)
#' X <- Prostate$X
#' y <- Prostate$y
#' cl <- makeCluster(4)
#' cvfit <- cv.ncvreg(X, y, cluster=cl, nfolds=length(y))}
#' 
#' # Survival
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#' 
#' cvfit <- cv.ncvsurv(X, y)
#' summary(cvfit)
#' plot(cvfit)
#' plot(cvfit, type="rsq")
#' @export cv.ncvreg

cv.ncvreg <- function(X, y, ..., cluster, nfolds=10, seed, fold, returnY=FALSE, trace=FALSE) {

  # Coercion
  if (!is.matrix(X)) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }

  fit <- ncvreg(X=X, y=y, ...)
  n <- length(y)
  E <- Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))
  if (fit$family=="binomial") {
    PE <- E
    if (!identical(sort(unique(y)), 0:1)) y <- as.double(y==max(y))
  }

  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  sde <- sqrt(.Machine$double.eps)
  if (missing(fold)) {
    if (fit$family=="binomial") {
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      fold1 <- 1:n1 %% nfolds
      fold0 <- (n1 + 1:n0) %% nfolds
      fold1[fold1==0] <- nfolds
      fold0[fold0==0] <- nfolds
      fold <- double(n)
      fold[y==1] <- sample(fold1)
      fold[y==0] <- sample(fold0)
    } else {
      fold <- sample(1:n %% nfolds)
      fold[fold==0] <- nfolds
    }
  } else {
    nfolds <- max(fold)
  }

  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$returnX <- FALSE
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("fold","fit","X", "y", "cv.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    fold.results <- parallel::parLapply(cl=cluster, X=1:nfolds, fun=cvf, XX=X, y=y, fold=fold, cv.args=cv.args)
  }

  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
    } else {
      if (trace) cat("Starting CV fold #", i, sep="","\n")
      res <- cvf(i, X, y, fold, cv.args)
    }
    E[fold==i, 1:res$nl] <- res$loss
    if (fit$family=="binomial") PE[fold==i, 1:res$nl] <- res$pe
    Y[fold==i, 1:res$nl] <- res$yhat
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, stats::sd) / sqrt(n)
  min <- which.min(cve)

  # Bias correction
  e <- sapply(1:nfolds, function(i) apply(E[fold==i, , drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  val <- list(cve=cve, cvse=cvse, fold=fold, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min],
              null.dev=mean(loss.ncvreg(y, rep(mean(y), n), fit$family)), Bias=Bias)
  if (fit$family=="binomial") {
    pe <- apply(PE, 2, mean)
    val$pe <- pe[is.finite(pe)]
  }
  if (returnY) val$Y <- Y
  structure(val, class="cv.ncvreg")
}
cvf <- function(i, XX, y, fold, cv.args) {
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  fit.i <- do.call("ncvreg", cv.args)

  X2 <- XX[fold==i, , drop=FALSE]
  y2 <- y[fold==i]
  yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
  loss <- loss.ncvreg(y2, yhat, fit.i$family)
  pe <- if (fit.i$family=="binomial") {(yhat < 0.5) == y2} else NULL
  list(loss=loss, pe=pe, nl=length(fit.i$lambda), yhat=yhat)
}
