summary.cv.ncvreg <- function(object, ...) {
  S <- pmax(object$null.dev - object$cve, 0)
  if (!('cv.ncvsurv' %in% class(object)) && object$fit$family=="gaussian") {
    rsq <- pmin(pmax(1 - object$cve/object$null.dev, 0), 1)
  } else {
    rsq <- pmin(pmax(1 - exp(object$cve-object$null.dev), 0), 1)
  }
  snr <- rsq/(1-rsq)
  nvars <- predict(object$fit, type="nvars")
  if ('cv.ncvsurv' %in% class(object)) {
    model <- 'Cox'
  } else {
    model <- switch(object$fit$family, gaussian="linear", binomial="logistic", poisson="Poisson")
  }
  val <- list(penalty=object$fit$penalty, model=model, n=object$fit$n, p=nrow(object$fit$beta)-1, min=object$min, lambda=object$lambda, cve=object$cve, r.squared=rsq, snr=snr, nvars=nvars)
  if (!('cv.ncvsurv' %in% class(object)) && object$fit$family=="gaussian") val$sigma <- sqrt(object$cve)
  if (!('cv.ncvsurv' %in% class(object)) && object$fit$family=="binomial") val$pe <- object$pe
  structure(val, class="summary.cv.ncvreg")
}
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
}
