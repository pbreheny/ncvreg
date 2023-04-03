#' Permutation fitting for ncvreg
#' 
#' Fits multiple penalized regression models in which the outcome is randomly
#' permuted, thereby allowing estimation of the marginal false discovery rate.
#' 
#' The function fits a penalized regression model to the actual data, then
#' repeats the process \code{N} times with a permuted version of the response
#' vector.  This allows estimation of the expected number of variables included
#' by chance for each value of \code{lambda}.  The ratio of this expected
#' quantity to the number of selected variables using the actual (non-permuted)
#' response is called the marginal false discovery rate (mFDR).
#' 
#' @param X The design matrix, without an intercept, as in \code{ncvreg}.
#' @param y The response vector, as in \code{ncvreg}.
#' @param ... Additional arguments to \code{ncvreg}.
#' @param permute What to permute.  If \code{'outcome'}, the response vector,
#' \code{y}, is permuted.  If \code{'residuals'}, the residuals are permuted.
#' This is only available for linear regression (i.e., for
#' \code{family='gaussian'}).  Note that permuting the residuals may take a
#' long time, as the residuals differ for each value of \code{lambda}, so
#' separate permutations are required at every value of \code{lambda}.  See
#' also \code{\link{permres}}.
#' @param N The number of permutation replications.  Default is 10.
#' @param seed You may set the seed of the random number generator in order to
#' obtain reproducible results.
#' @param trace If set to TRUE, perm.ncvreg will inform the user of its
#' progress by announcing the beginning of each permutation fit. Default is
#' FALSE.
#' @return An object with S3 class \code{"perm.ncvreg"} containing:
#' \item{EF}{The number of variables selected at each value of \code{lambda},
#' averaged over the permutation fits.} \item{S}{The actual number of selected
#' variables for the non-permuted data.} \item{mFDR}{The estimated marginal
#' false discovery rate (\code{EF/S}).} \item{fit}{The fitted \code{ncvreg}
#' object for the original (non-permuted) data.} \item{loss}{The loss/deviance
#' for each value of \code{lambda}, averaged over the permutation fits.  This
#' is an estimate of the explanatory power of the model under null conditions,
#' and can be used to adjust the loss of the fitted model in a manner akin to
#' the idea of an adjusted R-squared in classical regression.}
#' @author Patrick Breheny <patrick-breheny@@uiowa.edu>
#' @seealso \code{\link{ncvreg}}, \code{\link{plot.mfdr}}, \code{\link{mfdr}}
#' @examples
#' 
#' # Linear regression --------------------------------------------------
#' data(Prostate)
#' pmfit <- perm.ncvreg(Prostate$X, Prostate$y)
#' 
#' op <- par(mfcol=c(2,2))
#' plot(pmfit)
#' plot(pmfit, type="EF")
#' plot(pmfit$fit)
#' lam <- pmfit$fit$lambda
#' 
#' pmfit.r <- perm.ncvreg(Prostate$X, Prostate$y, permute='residuals')
#' plot(pmfit.r, col="red")              # Permuting residuals is
#' lines(lam, pmfit$mFDR, col="gray60")  # less conservative
#' par(op)
#' 
#' # Logistic regression ------------------------------------------------
#' data(Heart)
#' pmfit <- perm.ncvreg(Heart$X, Heart$y, family="binomial")
#' 
#' op <- par(mfcol=c(2,2))
#' plot(pmfit)
#' plot(pmfit, type="EF")
#' plot(pmfit$fit)
#' par(op)
#' 
#' @export perm.ncvreg
perm.ncvreg <- function(X, y, ..., permute=c("outcome", "residuals"), N=10, seed, trace=FALSE) {
  permute <- match.arg(permute)
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  fit <- ncvreg(X=X, y=y, returnX=TRUE, ...)
  if (fit$family != 'gaussian' & permute=='residuals') stop(paste0("Cannot permute residuals with family = ", fit$family), call.=FALSE)
  S <- predict(fit, type="nvars")

  if (permute=="outcome") {
    pfit <- fit.perm.ncvreg(fit, fit$y, fit$lambda, N, max(S), trace)
    S.perm <- pfit$S.perm
    L.perm <- pfit$L.perm
    EF <- pmin(apply(S.perm, 2, mean, na.rm=TRUE), S)
    mFDR <- EF/S
    mFDR[S==0] <- 0
    loss <- apply(L.perm, 2, mean)
  } else {
    n.l <- length(fit$lambda)
    EF <- mFDR <- loss <- double(n.l)
    for (i in 1:n.l) {
      pres <- permres(fit, fit$lambda[i], N=N, seed=seed, trace=trace)
      EF[i] <- pres$EF
      mFDR[i] <- pres$mFDR
      loss[i] <- pres$loss
      names(EF) <- names(mFDR) <- names(loss) <- fit$lambda
    }
  }

  fit <- structure(fit[1:11], class="ncvreg") ## Don't return X, y, etc.
  structure(list(EF=EF, S=S, mFDR=mFDR, fit=fit, loss=loss), class=c("perm.ncvreg", "mfdr"))
}
fit.perm.ncvreg <- function(fit, y, lam, N, maxdf, trace) {
  n.l <- length(lam)
  p <- ncol(fit$X)
  S.perm <- L.perm <- matrix(NA, N, n.l)
  for (i in 1:N) {
    if (trace) cat("Starting permutation fit #", i, sep="","\n")
    if (fit$family=="gaussian") {
      res <- .Call("cdfit_gaussian", fit$X, sample(y), fit$penalty, lam, 0.001, as.integer(1000), as.double(fit$gamma), fit$penalty.factor, fit$alpha, as.integer(maxdf), as.integer(TRUE))
      b <- matrix(res[[1]], p, n.l)
      ind <- is.na(res[[4]])
      L.perm[i,] <- res[[2]]
      L.perm[i, ind] <- NA
    } else {
      res <- .Call("cdfit_glm", fit$X, sample(y), fit$family, fit$penalty, lam, 0.001, as.integer(1000), as.double(fit$gamma), fit$penalty.factor, fit$alpha, as.integer(maxdf), as.integer(TRUE), as.integer(FALSE))
      b <- matrix(res[[2]], p, n.l)
      ind <- is.na(res[[5]])
      L.perm[i,] <- res[[3]]
      L.perm[i, ind] <- NA
    }
    S.perm[i,] <- apply(b[, drop=FALSE]!=0, 2, sum)
    S.perm[i, ind] <- NA
  }
  list(S.perm=S.perm, L.perm=L.perm)
}
