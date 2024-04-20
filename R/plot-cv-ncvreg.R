#' Plots the cross-validation curve from a cv.ncvreg object
#' 
#' Plots the cross-validation curve from a `cv.ncvreg` or `cv.ncvsurv` object,
#' along with standard error bars.
#' 
#' Error bars representing approximate 68% confidence intervals are plotted
#' along with the estimates across values of `lambda`. For `rsq` and `snr`
#' applied to models other than linear regression, the Cox-Snell R-squared is used.
#' 
#' @param x          A `cv.ncvreg` or `cv.ncvsurv` object.
#' @param log.l      Should horizontal axis be on the log scale?  Default is TRUE.
#' @param type       What to plot on the vertical axis:
#'   * `cve` plots the cross-validation error (deviance)
#'   * `rsq` plots an estimate of the fraction of the deviance explained by the model (R-squared)
#'   * `snr` plots an estimate of the signal-to-noise ratio
#'   * `scale` plots, for `family="gaussian"`, an estimate of the scale parameter (standard deviation)
#'   * `pred` plots, for `family="binomial"`, the estimated prediction error
#'   * `all` produces all of the above
#' @param selected   If `TRUE` (the default), places an axis on top of the plot
#'   denoting the number of variables in the model (i.e., that have a nonzero
#'   regression coefficient) at that value of `lambda`.
#' @param vertical.line   If `TRUE` (the default), draws a vertical line at the
#'   value where cross-validaton error is minimized.
#' @param col             Controls the color of the dots (CV estimates).
#' @param \dots           Other graphical parameters to [plot()]
#' 
#' @author Patrick Breheny
#' 
#' @seealso [ncvreg()], [cv.ncvreg()]
#' 
#' @references
#' Breheny P and Huang J. (2011) Coordinate descent algorithms for nonconvex
#' penalized regression, with applications to biological feature selection.
#' *Annals of Applied Statistics*, **5**: 232-253. \doi{10.1214/10-AOAS388}
#' 
#' @examples
#' # Linear regression --------------------------------------------------
#' data(Prostate)
#' cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
#' plot(cvfit)
#' op <- par(mfrow=c(2,2))
#' plot(cvfit, type="all")
#' par(op)
#' 
#' # Logistic regression ------------------------------------------------
#' data(Heart)
#' cvfit <- cv.ncvreg(Heart$X, Heart$y, family="binomial")
#' plot(cvfit)
#' op <- par(mfrow=c(2,2))
#' plot(cvfit, type="all")
#' par(op)
#' 
#' # Cox regression -----------------------------------------------------
#' data(Lung)
#' cvfit <- cv.ncvsurv(Lung$X, Lung$y)
#' op <- par(mfrow=c(1,2))
#' plot(cvfit)
#' plot(cvfit, type="rsq")
#' par(op)
#' @export

plot.cv.ncvreg <- function(x, log.l=TRUE, type=c("cve", "rsq", "scale", "snr", "pred", "all"), selected=TRUE, vertical.line=TRUE, col="red", ...) {
  type <- match.arg(type)
  if (type=="all") {
    plot(x, log.l=log.l, type="cve", selected=selected, ...)
    plot(x, log.l=log.l, type="rsq", selected=selected, ...)
    plot(x, log.l=log.l, type="snr", selected=selected, ...)
    if (length(x$fit$family)) {
      if (x$fit$family == "binomial") plot(x, log.l=log.l, type="pred", selected=selected, ...)
      if (x$fit$family == "gaussian") plot(x, log.l=log.l, type="scale", selected=selected, ...)
    }
    return(invisible(NULL))
  }
  l <- x$lambda
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)

  ## Calculate y
  L.cve <- x$cve - x$cvse
  U.cve <- x$cve + x$cvse
  if (type=="cve") {
    y <- x$cve
    L <- L.cve
    U <- U.cve
    ylab <- "Cross-validation error"
  } else if (type=="rsq" | type == "snr") {
    if (length(x$fit$family) && x$fit$family=='gaussian') {
      rsq <- pmin(pmax(1 - x$cve/x$null.dev, 0), 1)
      rsql <- pmin(pmax(1 - U.cve/x$null.dev, 0), 1)
      rsqu <- pmin(pmax(1 - L.cve/x$null.dev, 0), 1)
    } else {
      rsq <- pmin(pmax(1 - exp(x$cve-x$null.dev), 0), 1)
      rsql <- pmin(pmax(1 - exp(U.cve-x$null.dev), 0), 1)
      rsqu <- pmin(pmax(1 - exp(L.cve-x$null.dev), 0), 1)
    }
    if (type == "rsq") {
      y <- rsq
      L <- rsql
      U <- rsqu
      ylab <- ~R^2
    } else if(type=="snr") {
      y <- rsq/(1-rsq)
      L <- rsql/(1-rsql)
      U <- rsqu/(1-rsqu)
      ylab <- "Signal-to-noise ratio"
    }
  } else if (type=="scale") {
    if (x$fit$family == "binomial") stop("Scale parameter for binomial family fixed at 1", call.=FALSE)
    y <- sqrt(x$cve)
    L <- sqrt(L.cve)
    U <- sqrt(U.cve)
    ylab <- ~hat(sigma)
  } else if (type=="pred") {
    y <- x$pe
    n <- x$fit$n
    CI <- sapply(y, function(x) {binom.test(x*n, n, conf.level=0.68)$conf.int})
    L <- CI[1,]
    U <- CI[2,]
    ylab <- "Prediction error"
  }

  ind <- if (type=="pred") which(is.finite(l[1:length(x$pe)])) else which(is.finite(l[1:length(x$cve)]))
  ylim <- if (is.null(x$cvse)) range(y[ind]) else range(c(L[ind], U[ind]))
  aind <- intersect(ind, which((U-L)/diff(ylim) > 1e-3))
  plot.args = list(x=l[ind], y=y[ind], ylim=ylim, xlab=xlab, ylab=ylab, type="n", xlim=rev(range(l[ind])), las=1)
  new.args = list(...)
  if (length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  if (vertical.line) abline(v=l[x$min], lty=2, lwd=.5)
  suppressWarnings(arrows(x0=l[aind], x1=l[aind], y0=L[aind], y1=U[aind], code=3, angle=90, col="gray80", length=.05))
  points(l[ind], y[ind], col=col, pch=19, cex=.5)
  if (selected) {
    n.s <- predict(x$fit, lambda=x$lambda, type="nvars")
    axis(3, at=l, labels=n.s, tick=FALSE, line=-0.5)
    mtext("Variables selected", cex=0.8, line=1.5)
  }
}
