#' Model predictions based on a fitted `ncvsurv` object.
#' 
#' Similar to other predict methods, this function returns predictions from a
#' fitted `ncvsurv` object.
#' 
#' Estimation of baseline survival function conditional on the estimated values
#' of `beta` is carried out according to the method described in Chapter
#' 4.3 of Kalbfleish and Prentice.  In particular, it agrees exactly the
#' results returned by `survfit.coxph(..., type='kalbfleisch-prentice')`
#' in the `survival` package.
#' 
#' @aliases predict.ncvsurv coef.ncvsurv
#' 
#' @param object Fitted `"ncvsurv"` model object.
#' @param X Matrix of values at which predictions are to be made.  Not used for
#'   `type="coefficients"` or for some of the `type` settings in
#'   `predict`.
#' @param lambda Values of the regularization parameter `lambda` at which
#'   predictions are requested.  For values of `lambda` not in the sequence
#'   of fitted models, linear interpolation is used.
#' @param which Indices of the penalty parameter `lambda` at which predictions
#'   are required. By default, all indices are returned. If `lambda` is
#'   specified, this will override `which`.
#' @param type Type of prediction:
#'   * `link` returns the linear predictors
#'   * `response` gives the risk (i.e., exp(link))
#'   * `survival` returns the estimated survival function
#'   * `median` estimates median survival times
#'   The other options are all identical to their [ncvreg()] counterparts:
#'   * `coefficients` returns the coefficients
#'   * `vars` returns a list containing the indices and names of the nonzero variables at each value of `lambda`
#'   * `nvars` returns the number of nonzero coefficients at each value of `lambda`.
#' @param \dots Not used.
#' 
#' @returns The object returned depends on type.
#' 
#' @author Patrick Breheny <patrick-breheny@@uiowa.edu>
#' 
#' @seealso [ncvsurv()]
#' 
#' @references
#' * Breheny P and Huang J. (2011) Coordinate descent algorithms for nonconvex
#'   penalized regression, with applications to biological feature selection.
#'   *Annals of Applied Statistics*, **5**: 232-253. \doi{10.1214/10-AOAS388}
#'   
#' * Kalbfleish JD and Prentice RL (2002). *The Statistical Analysis of Failure
#'   Time Data*, 2nd edition. Wiley.

#' @examples
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#' 
#' fit <- ncvsurv(X,y)
#' coef(fit, lambda=0.05)
#' head(predict(fit, X, type="link", lambda=0.05))
#' head(predict(fit, X, type="response", lambda=0.05))
#' 
#' # Survival function
#' S <- predict(fit, X[1,], type="survival", lambda=0.05)[[1]]
#' S(100)
#' S <- predict(fit, X, type="survival", lambda=0.05)
#' plot(S, xlim=c(0,200))
#' 
#' # Medians
#' predict(fit, X[1,], type="median", lambda=0.05)
#' M <- predict(fit, X, type="median")
#' M[1:10, 1:10]
#' 
#' # Nonzero coefficients
#' predict(fit, type="vars", lambda = c(0.1, 0.01))
#' predict(fit, type="nvars", lambda = c(0.1, 0.01))
#' @export

predict.ncvsurv <- function(
    object,
    X,
    type=c("link", "response", "survival", "hazard", "median", "coefficients", "vars", "nvars"),
    lambda,
    which = 1:length(object$lambda), ...) {
  type <- match.arg(type)
  if (type %in% c("coefficients", "vars", "nvars")) {
    return(predict.ncvreg(
      object = object,
      X = X,
      type = type,
      lambda = lambda,
      which = which,
      ...
    ))
  }
  if (!missing(lambda)) {
    ind <- stats::approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    x <- ind %% 1
    beta <- (1 - x) * object$beta[, l, drop = FALSE] + x * object$beta[, r, drop = FALSE]
    colnames(beta) <- lam_names(lambda)
  } else {
    beta <- object$beta[, which, drop = FALSE]
  }

  if (missing(X)) {
    eta <- matrix(0, 1, 1)
    warning('Returning "baseline" prediction; supply X for more interesting prediction')
  } else {
    eta <- X %*% beta
  }
  if (type == 'link') return(drop(eta))
  if (type == 'response') return(drop(exp(eta)))

  if (!missing(lambda)) {
    W <- (1 - x) * exp(object$linear.predictors)[, l, drop = FALSE] +
      x * exp(object$linear.predictors)[, r, drop = FALSE]
  } else {
    W <- exp(object$linear.predictors)[, which, drop = FALSE]
  }

  if (type %in% c('survival', 'hazard') & ncol(W) > 1) {
    stop('Can only return type="survival" for a single lambda value', call. = FALSE)
  } 
  if (type %in% c('survival', 'hazard')) val <- vector('list', length(eta))
  if (type == 'median') val <- matrix(NA, nrow(eta), ncol(eta))
  for (j in 1:ncol(eta)) {
    # Estimate baseline hazard
    w <- W[,j]
    r <- rev(cumsum(rev(w)))
    a <- ifelse(object$fail, (1 - w / r)^(1 / w), 1)
    S0 <- c(1, cumprod(a))
    H0 <- c(0, cumsum(1 - a))
    x <- c(0, object$time)
    for (i in 1:nrow(eta)) {
      S <- S0^exp(eta[i,j])
      H <- H0*exp(eta[i,j])
      if (type == 'survival') val[[i]] <- approxfun(x, S, method = 'constant', ties = "ordered")
      else if (type == 'hazard') val[[i]] <- approxfun(x, H, method = 'constant', ties = "ordered")
      else if (type == 'median') {
        if (any(S < 0.5)) {
          val[i,j] <- x[min(which(S < .5))]
        }
      }
    }
  }
  if (type %in% c('survival', 'hazard')) {
    class(val) <- c('ncvsurv.func', class(val))
    attr(val, 'time') <- object$time
  }
  if (type == 'median') val <- drop(val)
  val
}
