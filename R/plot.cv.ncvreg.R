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
  } else if (type=="snr") {
    S <- pmax(x$null.dev - x$cve, 0)
    y <- S/(x$cve)
    L <- S/U.cve
    U <- S/L.cve

  } else if (type=="scale") {
    if (x$fit$family == "binomial") stop("Scale parameter for binomial family fixed at 1")
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
  if (vertical.line) abline(v=l[x$min],lty=2,lwd=.5)
  suppressWarnings(arrows(x0=l[aind], x1=l[aind], y0=L[aind], y1=U[aind], code=3, angle=90, col="gray80", length=.05))
  points(l[ind], y[ind], col=col, pch=19, cex=.5)
  if (selected) {
    n.s <- predict(x$fit, lambda=x$lambda, type="nvars")
    axis(3, at=l, labels=n.s, tick=FALSE, line=-0.5)
    mtext("Variables selected", cex=0.8, line=1.5)
  }
}
