#' Model predictions based on a fitted ncvreg object.
#' 
#' Similar to other predict methods, this function returns predictions from a
#' fitted \code{ncvreg} object.
#' 
#' 
#' @aliases predict.ncvreg coef.ncvreg
#' @param object Fitted \code{ncvreg} model object.
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
#' \code{"response"} gives the fitted values; \code{"class"} returns the
#' binomial outcome with the highest probability; \code{"coefficients"} returns
#' the coefficients; \code{"vars"} returns a list containing the indices and
#' names of the nonzero variables at each value of \code{lambda};
#' \code{"nvars"} returns the number of nonzero coefficients at each value of
#' \code{lambda}.
#' @param drop If coefficients for a single value of \code{lambda} are to be
#' returned, reduce dimensions to a vector?  Setting \code{drop=FALSE} returns
#' a 1-column matrix.
#' @param \dots Not used.
#' @return The object returned depends on type.
#' @author Patrick Breheny
#' @seealso \code{\link{ncvreg}}
#' @references Breheny P and Huang J. (2011) Coordinate descentalgorithms for
#' nonconvex penalized regression, with applications to biological feature
#' selection.  \emph{Annals of Applied Statistics}, \strong{5}: 232-253.
#' c("\\Sexpr[results=rd]{tools:::Rd_expr_doi(\"#1\")}",
#' "10.1214/10-AOAS388")\Sexpr{tools:::Rd_expr_doi("10.1214/10-AOAS388")}
#' @examples
#' 
#' data(Heart)
#' 
#' fit <- ncvreg(Heart$X, Heart$y, family="binomial")
#' coef(fit, lambda=0.05)
#' head(predict(fit, Heart$X, type="link", lambda=0.05))
#' head(predict(fit, Heart$X, type="response", lambda=0.05))
#' head(predict(fit, Heart$X, type="class", lambda=0.05))
#' predict(fit, type="vars", lambda=c(0.05, 0.01))
#' predict(fit, type="nvars", lambda=c(0.05, 0.01))
#' @export

predict.ncvreg <- function(object, X, type=c("link", "response", "class", "coefficients", "vars", "nvars"),
                           lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)
  beta <- coef.ncvreg(object, lambda=lambda, which=which, drop=FALSE)
  if (type=="coefficients") return(beta)
  if (!inherits(object, 'ncvsurv')) {
    alpha <- beta[1,]
    beta <- beta[-1, , drop=FALSE]
  } else {
    beta <- beta
  }

  if (type=="nvars") return(apply(beta!=0, 2, sum))
  if (type=="vars") return(drop(apply(beta!=0, 2, FUN=which)))
  eta <- sweep(X %*% beta, 2, alpha, "+")
  if (type=="link" || object$family=="gaussian") return(drop(eta))
  resp <- switch(object$family,
                 binomial = exp(eta)/(1+exp(eta)),
                 poisson = exp(eta))
  if (type=="response") return(drop(resp))
  if (type=="class") {
    if (object$family=="binomial") {
      return(drop(1*(eta>0)))
    } else {
      stop("type='class' can only be used with family='binomial'", call.=FALSE)
    }
  }
}

#' @rdname predict.ncvreg
#' @export
coef.ncvreg <- function(object, lambda, which=1:length(object$lambda), drop=TRUE, ...) {
  if (!missing(lambda)) {
    if (max(lambda) > max(object$lambda) | min(lambda) < min(object$lambda)) {
      stop('Supplied lambda value(s) are outside the range of the model fit.', call.=FALSE)
    }
    ind <- stats::approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    beta <- (1-w)*object$beta[, l, drop=FALSE] + w*object$beta[, r, drop=FALSE]
    colnames(beta) <- lam_names(lambda)
  }
  else beta <- object$beta[, which, drop=FALSE]
  if (drop) return(drop(beta)) else return(beta)
}
