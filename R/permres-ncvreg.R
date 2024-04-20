#' Permute residuals for a fitted ncvreg model
#' 
#' Fits multiple penalized regression models in which the residuals are
#' randomly permuted, thereby allowing estimation of the marginal false
#' discovery rate.
#' 
#' The function fits a penalized regression model to the actual data, then
#' repeats the process `N` times with a permuted version of the response
#' vector. This allows estimation of the expected number of variables included
#' by chance for each value of `lambda`. The ratio of this expected
#' quantity to the number of selected variables using the actual (non-permuted)
#' response is called the marginal false discovery rate (mFDR).
#' 
#' @aliases permres permres.ncvreg
#' 
#' @param fit A fitted ncvreg model, as produced by [ncvreg()]. To use with
#'   `permres`, the model must be fit using the `returnX=TRUE` option.
#' @param lambda The regularization parameter to use for estimating residuals.
#'   Unlike [perm.ncvreg()], `permres()` calculates EF and mFDR for a specific 
#'   `lambda` value, not an entire path. As a result, it runs much faster.
#' @param N The number of permutation replications.  Default is 10.
#' @param seed You may set the seed of the random number generator in order to
#'   obtain reproducible results.
#' @param trace If set to TRUE, perm.ncvreg will inform the user of its
#'   progress by announcing the beginning of each permutation fit. Default is
#'   FALSE.
#' @param \dots Not used.
#' 
#' @returns A list with the following components:
#' \item{EF}{The number of variables selected at each value of `lambda`, averaged over the permutation fits.}
#' \item{S}{The actual number of selected variables for the non-permuted data.}
#' \item{mFDR}{The estimated marginal false discovery rate (`EF/S`).}
#' \item{loss}{The loss/deviance, averaged over the permutation fits. This is an estimate of the explanatory power of the model under null conditions, and can be used to adjust the loss of the fitted model in a manner akin to the idea of an adjusted R-squared in classical regression.}
#' 
#' @author Patrick Breheny <patrick-breheny@@uiowa.edu>
#' 
#' @seealso [ncvreg()], `[mfdr()], [perm.ncvreg()]
#' 
#' @examples
#' data(Prostate)
#' fit <- ncvreg(Prostate$X, Prostate$y, N=50)
#' permres(fit, lambda=0.15)
#' @export permres

permres <- function(fit, ...) UseMethod("permres")

#' @rdname permres
#' @export
permres.ncvreg <- function(fit, lambda, N=10, seed, trace=FALSE, ...) {
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  if (!is.numeric(lambda)) stop("lambda must be numeric", call.=FALSE)
  if (length(lambda) != 1) stop("lambda must be a single number, not a vector", call.=FALSE)
  if (is.null(fit$X)) stop("must run ncvreg with returnX=TRUE", call.=FALSE)

  S <- predict(fit, type="nvars", lambda=lambda)
  y <- fit$y - predict(fit, fit$X, lambda=lambda)
  y <- y - mean(y)
  lam <- c(fit$lambda[fit$lambda > lambda], lambda)
  pfit <- fit.perm.ncvreg(fit, y, lam, N, S, trace)
  S.perm <- pfit$S.perm[, ncol(pfit$S.perm)]
  L.perm <- pfit$L.perm[, ncol(pfit$S.perm)]

  EF <- mean(S.perm, na.rm=TRUE)
  mFDR <- EF/S
  mFDR[S==0] <- 0
  list(EF=EF, S=S, mFDR=mFDR, loss=mean(L.perm, na.rm=TRUE))
}
