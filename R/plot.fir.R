plot.fir <- function(x, type=c("FIR", "EF"), log.l=FALSE, legend=TRUE, ...) {
  type <- match.arg(type)
  l <- as.numeric(names(x$FIR))
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)
  
  if (type=="FIR") {
    plot.args <- list(x=l, y=x$FIR, xlim=rev(range(l)), las=1, xlab=xlab, ylab="FIR", type="l")
    new.args <- list(...)
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
  } else if (type=="EF") {
    col <- c("#FF4E37FF", "#008DFFFF")
    plot.args <- list(x=l, y=cbind(x$S, x$EF), xlim=rev(range(l)), las=1, lty=1, lwd=2, xlab=xlab, ylab="# of variables", type="l", col=col)
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
transform.coord <- function(x,p) {
  ba <- (x[2]-x[1])/(p[2]-p[1])
  a <- x[1]-p[1]*ba
  b <- a + ba
  c(a,b)
}
