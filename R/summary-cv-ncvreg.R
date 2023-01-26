#' Summarizing cross-validation-based inference
#' 
#' Summary method for \code{cv.ncvreg} objects
#' 
#' 
#' @aliases summary.cv.ncvreg print.summary.cv.ncvreg
#' @param object A \code{"cv.ncvreg"} or \code{"cv.ncvsurv"} object.
#' @param x A \code{"summary.cv.ncvreg"} object.
#' @param digits Number of digits past the decimal point to print out.  Can be
#' a vector specifying different display digits for each of the five
#' non-integer printed values.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{summary.cv.ncvreg} produces an object with S3 class
#' \code{"summary.cv.ncvreg"}.  The class has its own print method and contains
#' the following list elements: \describe{ \item{penalty}{The penalty used by
#' \code{ncvreg}.} \item{model}{Either \code{"linear"} or \code{"logistic"},
#' depending on the \code{family} option in \code{ncvreg}.} \item{n}{Number of
#' observations} \item{p}{Number of regression coefficients (not including the
#' intercept).} \item{min}{The index of \code{lambda} with the smallest
#' cross-validation error.} \item{lambda}{The sequence of \code{lambda} values
#' used by \code{cv.ncvreg}.} \item{cve}{Cross-validation error (deviance).}
#' \item{r.squared}{Proportion of variance explained by the model, as estimated
#' by cross-validation.  For models outside of linear regression, the Cox-Snell
#' approach to defining R-squared is used.} \item{snr}{Signal to noise ratio,
#' as estimated by cross-validation.} \item{sigma}{For linear regression
#' models, the scale parameter estimate.} \item{pe}{For logistic regression
#' models, the prediction error (misclassification error).} }
#' @author Patrick Breheny
#' @seealso \code{\link{ncvreg}}, \code{\link{cv.ncvreg}},
#' \code{\link{plot.cv.ncvreg}}
#' @references Breheny P and Huang J. (2011) Coordinate descentalgorithms for
#' nonconvex penalized regression, with applications to biological feature
#' selection.  \emph{Annals of Applied Statistics}, \strong{5}: 232-253.
#' c("\\Sexpr[results=rd]{tools:::Rd_expr_doi(\"#1\")}",
#' "10.1214/10-AOAS388")\Sexpr{tools:::Rd_expr_doi("10.1214/10-AOAS388")}
#' @examples
#' 
#' # Linear regression --------------------------------------------------
#' data(Prostate)
#' cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
#' summary(cvfit)
#' 
#' # Logistic regression ------------------------------------------------
#' data(Heart)
#' cvfit <- cv.ncvreg(Heart$X, Heart$y, family="binomial")
#' summary(cvfit)
#' 
#' # Cox regression -----------------------------------------------------
#' data(Lung)
#' cvfit <- cv.ncvsurv(Lung$X, Lung$y)
#' summary(cvfit)
#' @export

summary.cv.ncvreg <- function(object, ...) {
  S <- pmax(object$null.dev - object$cve, 0)
  if (!inherits(object, 'cv.ncvsurv') && object$fit$family=="gaussian") {
    rsq <- pmin(pmax(1 - object$cve/object$null.dev, 0), 1)
  } else {
    rsq <- pmin(pmax(1 - exp(object$cve-object$null.dev), 0), 1)
  }
  snr <- rsq/(1-rsq)
  nvars <- predict(object$fit, type="nvars")
  if (inherits(object, 'cv.ncvsurv')) {
    model <- 'Cox'
  } else {
    model <- switch(object$fit$family, gaussian="linear", binomial="logistic", poisson="Poisson")
  }
  val <- list(penalty=object$fit$penalty, model=model, n=object$fit$n, p=nrow(object$fit$beta)-1, 
              min=object$min, lambda=object$lambda, cve=object$cve, r.squared=rsq, snr=snr, 
              nvars=nvars)
  if (!inherits(object, 'cv.ncvsurv') && object$fit$family=="gaussian") {
    val$sigma <- sqrt(object$cve)
    val$fit_summary <- summary(object$fit, object$lambda.min, sigma=val$sigma[object$min], ...)
  } else {
    val$fit_summary <- summary(object$fit, object$lambda.min, ...)
  }
  if (!inherits(object, 'cv.ncvsurv') && object$fit$family=="binomial") val$pe <- object$pe
  structure(val, class="summary.cv.ncvreg")
}

#' @rdname summary.cv.ncvreg
#' @export
print.summary.cv.ncvreg <- function(x, digits, ...) {
  digits <- if (missing(digits)) digits <- c(2, 4, 2, 2, 3) else rep(digits, length.out=5)
  cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, ", p=", x$p, "\n", sep="")
  cat("At minimum cross-validation error (lambda=", formatC(x$lambda[x$min], digits[2], format="f"), "):\n", sep="")
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", x$nvars[x$min], "\n", sep="")
  cat("  Cross-validation error (deviance): ", formatC(min(x$cve), digits[1], format="f"), "\n", sep="")
  cat("  R-squared: ", formatC(max(x$r.squared), digits[3], format="f"), "\n", sep="")
  cat("  Signal-to-noise ratio: ", formatC(max(x$snr), digits[4], format="f"), "\n", sep="")
  if (x$model == "logistic") cat("  Prediction error: ", formatC(x$pe[x$min], digits[5], format="f"), "\n", sep="")
  if (x$model == "linear") cat("  Scale estimate (sigma): ", formatC(sqrt(x$cve[x$min]), digits[5], format="f"), "\n", sep="")
  if (x$nvars[x$min] > 0) print.summary.ncvreg(x$fit_summary, ...)
}
