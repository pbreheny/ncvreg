#' Plot marginal false discovery rate curves
#' 
#' Plot marginal false discovery rate curves from an \code{"mfdr"} or
#' \code{"perm.ncvreg"} object.
#' 
#' 
#' @param x A \code{"perm.ncvreg"} or \code{"mfdr"} object.
#' @param type What to plot on the vertical axis.  \code{mFDR} plots the
#' marginal false discovery rate; \code{EF} plots the expected number of false
#' discoveries along with the actual number of variables included in the model.
#' @param log.l Should horizontal axis be on the log scale?  Default is FALSE.
#' @param selected If \code{TRUE} (the default), places an axis on top of the
#' plot denoting the number of variables in the model (i.e., that have a
#' nonzero regression coefficient) at that value of \code{lambda}.
#' @param legend For \code{type="EF"} plots, draw a legend to indicate which
#' line is for the actual selections and which line is for the expected number
#' of false discoveries?  Default is \code{TRUE}.
#' @param \dots Other graphical parameters to pass to \code{plot}
#' @author Patrick Breheny
#' @seealso \code{\link{mfdr}}, \code{\link{perm.ncvreg}}
#' @references Breheny P (2019). Marginal false discovery rates for penalized
#' regression models. Biostatistics, 20: 299-314.
#' @examples
#' 
#' data(Prostate)
#' fit <- ncvreg(Prostate$X, Prostate$y)
#' 
#' obj <- mfdr(fit)
#' obj[1:10,]
#' 
#' # Some plotting options
#' plot(obj)
#' plot(obj, type="EF")
#' plot(obj, log=TRUE)
#' 
#' 
#' # Comparison with perm.ncvreg
#' op <- par(mfrow=c(2,2))
#' plot(obj)
#' plot(obj, type="EF")
#' pmfit <- perm.ncvreg(Prostate$X, Prostate$y)
#' plot(pmfit)
#' plot(pmfit, type="EF")
#' par(op)
#' @export

plot.mfdr <- function(x, type=c("mFDR", "EF"), log.l=FALSE, selected=TRUE, legend=TRUE, ...) {
  if (inherits(x, "perm.ncvreg")) {
    l <- x$fit$lambda
    x <- data.frame(EF=x$EF, S=x$S, mFDR=x$mFDR)
  } else {
    l <- as.double(rownames(x))
  }
  type <- match.arg(type)
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)

  if (type=="mFDR") {
    plot.args <- list(x=l, y=x$mFDR, xlim=rev(range(l)), las=1, xlab=xlab, ylab="mFDR", type="l")
    new.args <- list(...)
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
    if (selected) {
      ll <- seq(min(l), max(l), len=6)
      ind <- sapply(ll, function(x) which.min(abs(x-l)))
      axis(3, at = ll, labels = x$S[ind], tick = TRUE)
      mtext("Variables selected", cex = 0.8, line = 2)
    }
  } else if (type=="EF") {
    col <- c("#FF4E37FF", "#008DFFFF")
    plot.args <- list(x=l, y=cbind(x$S, x$EF), xlim=rev(range(l)), las=1, lty=1, lwd=2, xlab=xlab, ylab="# of variables", type=c("S", "l"), col=col)
    new.args <- list(...)
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("matplot", plot.args)
    if (legend) {
      x <- mean(par("usr")[1:2])
      yy <- transform.coord(par("usr")[3:4], par("plt")[3:4])
      y <- mean(c(yy[2], par("usr")[4]))
      legend(legend=c("Selected", "E(False)"), col=col, x, y, xpd = NA, bty = "n", xjust = 0.5, yjust = 0.5, horiz = TRUE, lwd=plot.args$lwd)
    }
  }
}
transform.coord <- function(x, p) {
  ba <- (x[2]-x[1])/(p[2]-p[1])
  a <- x[1]-p[1]*ba
  b <- a + ba
  c(a, b)
}
