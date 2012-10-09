plot.cv.ncvreg <- function(x, log.l=TRUE, ...)
{
  L <- x$cve - x$cvse
  U <- x$cve + x$cvse
  l <- x$lambda
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)
  
  ind <- which(((U - L) > 1e-4) & is.finite(l))
  plot.args = list(x=l[ind], y=x$cve[ind], ylim=range(c(L[ind],U[ind])), xlab=xlab, ylab="Cross-validation error", type="n", xlim=rev(range(l[ind])))
  new.args = list(...)
  if (length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  abline(v=l[x$min],lty=2,lwd=.5)
  arrows(x0=l[ind], x1=l[ind], y0=L[ind], y1=U[ind], code=3, angle=90, col="gray80", length=.05)
  points(l[ind], x$cve[ind], col="red", pch=19, cex=.5)
}
