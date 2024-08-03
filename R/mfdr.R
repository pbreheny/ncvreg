#' Marginal false discovery rates
#' 
#' Estimates the marginal false discovery rate (mFDR) of a penalized regression
#' model.
#' 
#' The function estimates the marginal false discovery rate (mFDR) for a
#' penalized regression model.  The estimate tends to be accurate in most
#' settings, but will be slightly conservative if predictors are highly
#' correlated.  For an alternative way of estimating the mFDR, typically more
#' accurate in highly correlated cases, see [perm.ncvreg()].
#' 
#' @param fit   An `ncvreg` or `ncvsurv` object.
#' @param X     The model matrix corresponding to `fit`. This is not necessary
#'   for linear regression, but in logistic and Cox regression, the mFDR depends
#'   on X. It is not necessary to supply `X` if it is already contained in `fit`;
#'   i.e., if `ncvreg`/`ncvsurv` was run with `returnX=TRUE`.
#' 
#' @returns An object with S3 class `mfdr` inheriting from `data.frame`, containing:
#' \describe{
#'   \item{EF}{The number of variables selected at each value of `lambda`, averaged over the permutation fits.}
#'   \item{S}{The actual number of selected variables for the non-permuted data.}
#'   \item{mFDR}{The estimated marginal false discovery rate (`EF/S`).}
#' }
#' 
#' @author Patrick Breheny and Ryan Miller
#' 
#' @seealso [ncvreg()], [ncvsurv()], [plot.mfdr()], [perm.ncvreg()]
#' 
#' @examples
#' # Linear regression --------------------------------
#' data(Prostate)
#' fit <- ncvreg(Prostate$X, Prostate$y)
#' 
#' obj <- mfdr(fit)
#' obj[1:10,]
#' 
#' # Comparison with perm.ncvreg
#' op <- par(mfrow=c(2,2))
#' plot(obj)
#' plot(obj, type="EF")
#' pmfit <- perm.ncvreg(Prostate$X, Prostate$y)
#' plot(pmfit)
#' plot(pmfit, type="EF")
#' par(op)
#' 
#' # Logistic regression ------------------------------
#' data(Heart)
#' fit <- ncvreg(Heart$X, Heart$y, family="binomial")
#' obj <- mfdr(fit)
#' head(obj)
#' op <- par(mfrow=c(1,2))
#' plot(obj)
#' plot(obj, type="EF")
#' par(op)
#' 
#' # Cox regression -----------------------------------
#' data(Lung)
#' fit <- ncvsurv(Lung$X, Lung$y)
#' obj <- mfdr(fit)
#' head(obj)
#' op <- par(mfrow=c(1,2))
#' plot(obj)
#' plot(obj, type="EF")
#' par(op)
#' @export mfdr

mfdr <- function(fit, X) {
  # Setup
  if (!inherits(fit, 'ncvreg')) stop('"fit" must be an ncvreg object', call.=FALSE)
  linear <- !inherits(fit, "ncvsurv") && fit$family == "gaussian"
  S0 <- sum(fit$penalty.factor==0)
  S <- predict(fit, type="nvars") - S0
  if (inherits(fit, "ncvsurv") || fit$family == "binomial") {
    if (!("X" %in% names(fit))) {
      if (missing(X)) {
        stop("For Cox/GLM models, you must either supply X or run ncvsurv with returnX=TRUE to calculate an mFDR", call.=FALSE)
      } else {
        if (inherits(fit, "ncvsurv")) {
          fit$X <- std(X)[fit$order,]
        } else {
          fit$X <- std(X)
        }
      }
    }
  }

  # Call C functions
  if (inherits(fit, "ncvsurv")) {
    EF <- .Call("mfdr_cox", fit)
  } else {
    if (fit$family == "binomial") {
      EF <- .Call("mfdr_binomial", fit)
    } else if (fit$family == "gaussian") {
      EF <- .Call("mfdr_gaussian", fit)
    }
  }

  # Calculate rate, return
  EF <- pmin(EF - S0, S)
  mFDR <- EF/S
  mFDR[S==0] <- 0
  df <- data.frame(EF=EF, S=S, mFDR=mFDR)
  rownames(df) <- lam_names(fit$lambda)
  structure(df, class=c("mfdr", "data.frame"))
}
