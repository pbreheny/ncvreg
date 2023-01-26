#' @keywords internal 
#' @aliases ncvreg-package NULL
#' @references Breheny P and Huang J. (2011) Coordinate descent algorithms for
#' nonconvex penalized regression, with applications to biological feature
#' selection. *Annals of Applied Statistics*, **5**: 232-253.
#' @examples
#' \donttest{vignette("getting-started", package="ncvreg")}
"_PACKAGE"

#' @useDynLib ncvreg, .registration = TRUE
#' @import stats
#' @import graphics
NULL

#' Internal ncvreg functions
#' 
#' Internal ncvreg functions
#' 
#' These are not intended for use by users. \code{convexMin} calculates the
#' lowest index for which the penalized objective function is locally convex.
#' \code{setupLambda} creates an appropriate vector of regularization parameter
#' values.  \code{loss.ncvreg} calculates the value of the loss function for
#' the given predictions (used for cross-validation).
#' 
#' @aliases setupLambda convexMin
#' @author Patrick Breheny <patrick-breheny@@uiowa.edu>
#' @keywords internal
