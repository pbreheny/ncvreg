#' Model predictions based on a fitted "ncvsurv" object.
#' 
#' Similar to other predict methods, this function returns predictions from a
#' fitted \code{"ncvsurv"} object.
#' 
#' Estimation of baseline survival function conditional on the estimated values
#' of \code{beta} is carried out according to the method described in Chapter
#' 4.3 of Kalbfleish and Prentice.  In particular, it agrees exactly the
#' results returned by \code{survfit.coxph(..., type='kalbfleisch-prentice')}
#' in the \code{survival} package.
#' 
#' @aliases predict.ncvsurv coef.ncvsurv
#' @param object Fitted \code{"ncvsurv"} model object.
#' @param X Matrix of values at which predictions are to be made.  Not used for
#' \code{type="coefficients"} or for some of the \code{type} settings in
#' \code{predict}.
#' @param lambda Values of the regularization parameter \code{lambda} at which
#' predictions are requested.  For values of \code{lambda} not in the sequence
#' of fitted models, linear interpolation is used.
#' @param which Indices of the penalty parameter \code{lambda} at which
#' predictions are required.  By default, all indices are returned.  If
#' \code{lambda} is specified, this will override \code{which}.
#' @param type Type of prediction: \code{"link"} returns the linear predictors;
#' \code{"response"} gives the risk (i.e., exp(link)); \code{"survival"}
#' returns the estimated survival function; \code{"median"} estimates median
#' survival times.  The other options are all identical to their \code{ncvreg}
#' counterparts: \code{"coefficients"} returns the coefficients; \code{"vars"}
#' returns a list containing the indices and names of the nonzero variables at
#' each value of \code{lambda}; \code{"nvars"} returns the number of nonzero
#' coefficients at each value of \code{lambda}.
#' @param \dots Not used.
#' @return The object returned depends on type.
#' @author Patrick Breheny <patrick-breheny@@uiowa.edu>
#' @seealso \code{\link{ncvsurv}}
#' @references \itemize{ \item Breheny P and Huang J. (2011) Coordinate
#' descentalgorithms for nonconvex penalized regression, with applications to
#' biological feature selection.  \emph{Annals of Applied Statistics},
#' \strong{5}: 232-253.  c("\\Sexpr[results=rd]{tools:::Rd_expr_doi(\"#1\")}",
#' "10.1214/10-AOAS388")\Sexpr{tools:::Rd_expr_doi("10.1214/10-AOAS388")}
#' 
#' \item Kalbfleish JD and Prentice RL (2002). \emph{The Statistical Analysis
#' of Failure Time Data}, 2nd edition. Wiley.  }
#' @examples
#' 
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
#' S <- predict(fit, X[1,], type="survival", lambda=0.05)
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
#' predict(fit, type="vars", lambda=c(0.1, 0.01))
#' predict(fit, type="nvars", lambda=c(0.1, 0.01))
#' @export

predict.ncvsurv <- function(object, X, type=c("link", "response", "survival",
                                              "median", "coefficients", "vars",
                                              "nvars"),
                            lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)
  if (type %in% c("coefficients", "vars", "nvars")) {
    return(predict.ncvreg(object=object, X=X, type=type, lambda=lambda, which=which, ...))
  }
  if (!missing(lambda)) {
    ind <- stats::approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    x <- ind %% 1
    beta <- (1-x)*object$beta[, l, drop=FALSE] + x*object$beta[, r, drop=FALSE]
    colnames(beta) <- lamNames(lambda)
  } else {
    beta <- object$beta[, which, drop=FALSE]
  }

  eta <- X %*% beta
  if (type=='link') return(drop(eta))
  if (type=='response') return(drop(exp(eta)))

  if (!missing(lambda)) {
    W <- (1-x)*exp(object$linear.predictors)[, l, drop=FALSE] + x*exp(object$linear.predictors)[, r, drop=FALSE]
  } else {
    W <- exp(object$linear.predictors)[, which, drop=FALSE]
  }
  if (type == 'survival' & ncol(W) > 1) stop('Can only return type="survival" for a single lambda value', call.=FALSE)
  if (type == 'survival') val <- vector('list', length(eta))
  if (type == 'median') val <- matrix(NA, nrow(eta), ncol(eta))
  for (j in 1:ncol(eta)) {
    # Estimate baseline hazard
    w <- W[,j]
    r <- rev(cumsum(rev(w)))
    a <- ifelse(object$fail, (1-w/r)^(1/w), 1)
    S0 <- c(1, cumprod(a))
    x <- c(0, object$time)
    for (i in 1:nrow(eta)) {
      S <- S0^exp(eta[i,j])
      if (type == 'survival') val[[i]] <- stats::approxfun(x, S, method='constant', ties=function(x) utils::tail(x, 1))
      if (type == 'median') {
        if (any(S < 0.5)) {
          val[i,j] <- x[min(which(S < .5))]
        }
      }
    }
  }
  if (type == 'survival') {
    if (nrow(eta)==1) val <- val[[1]]
    class(val) <- c('ncvsurv.func', class(val))
    attr(val, 'time') <- object$time
  }
  if (type == 'median') val <- drop(val)
  val
}
