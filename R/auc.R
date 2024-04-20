#' AUC for cv.ncvsurv objects
#' 
#' Calculates the cross-validated AUC (concordance) from a `cv.ncvsurv` object.
#' 
#' The area under the curve (AUC), or equivalently, the concordance statistic
#' (C), is calculated according to the procedure described in van Houwelingen
#' and Putter (2011). The function calls [survival::concordancefit()], except
#' cross-validated linear predictors are used to guard against overfitting.
#' Thus, the values returned by `AUC.cv.ncvsurv()` will be lower than those you
#' would obtain with `concordancefit()` if you fit the full (unpenalized) model.
#' 
#' @param obj   A `cv.ncvsurv` object. You must run [cv.ncvsurv()] with the
#'   option `returnY=TRUE` in order for `AUC()` to work.
#' @param ...   For S3 method compatibility; not used
#' 
#' @aliases AUC 
#' 
#' @references van Houwelingen H, Putter H (2011). Dynamic Prediction in
#' Clinical Survival Analysis.  CRC Press.
#' 
#' @author Patrick Breheny, Brandon Butcher, and Lawrence Hunsicker
#' 
#' @seealso [cv.ncvsurv()], [survival::concordancefit()]
#' 
#' @rdname AUC
#' 
#' @examples
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#' 
#' cvfit <- cv.ncvsurv(X, y, returnY=TRUE)
#' head(AUC(cvfit))
#' lam <- cvfit$lambda
#' plot(lam, AUC(cvfit), xlim=rev(range(lam)), lwd=3, type='l',
#'      las=1, xlab=expression(lambda), ylab='AUC')
#' @export

AUC.cv.ncvsurv <- function(obj, ...) {
  if (!("Y" %in% names(obj))) stop("Must run cv.ncvsurv with 'returnY=TRUE' in order to calculate AUC", call.=FALSE)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is needed for AUC() to work. Please install it.", call. = FALSE)
  }  
  if (utils::packageVersion("survival") < "3.2.10") stop("AUC.cv.ncvsurv requires version 3.2.10 of 'survival' package or higher", call.=FALSE)
  S <- survival::Surv(obj$fit$time, obj$fit$fail)
  apply(obj$Y[obj$fit$order, ], 2, concord, y = S)
}

#' @export
AUC <- function(obj, ...) UseMethod("AUC")
concord <- function(x, y) {
  survival::concordancefit(y, -x)$concordance
}
