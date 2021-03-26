AUC.cv.ncvsurv <- function(obj, ...) {
  if (!("Y" %in% names(obj))) stop("Must run cv.ncvsurv with 'returnY=TRUE' in order to calculate AUC", call.=FALSE)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is needed for AUC() to work. Please install it.", call. = FALSE)
  }  
  if (packageVersion("survival") < "3.2.10") stop("AUC.cv.ncvsurv requires version 3.2.10 of 'survival' package or higher", call.=FALSE)
  S <- survival::Surv(obj$fit$time, obj$fit$fail)
  res <- apply(obj$Y, 2, concord, y = S)
  num <- res['concordant',] + 0.5*res['tied.x',] + 0.5*res['tied.y',] + 0.5*res['tied.xy',]
  num/sum(res[,1])
}
AUC <- function(obj, ...) UseMethod("AUC")
concord <- function(x, y) {
  survival::concordancefit(y, -x)$count
}
