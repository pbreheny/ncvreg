AUC.cv.ncvsurv <- function(obj, ...) {
  if (!("Y" %in% names(obj))) stop("Must run cv.ncvsurv with 'returnY=TRUE' in order to calculate AUC")
  SURV <- get("Surv", asNamespace("survival"))
  S <- SURV(obj$fit$time, obj$fit$fail)
  CONC <- get("survConcordance.fit", asNamespace("survival"))
  res <- apply(obj$Y, 2, CONC, y = S)
  num <- res['concordant',] + 0.5*res['tied.risk',] + 0.5*res['tied.time',]
  num/sum(res[1:4,1])
}
AUC <- function(obj, ...) UseMethod("AUC")
