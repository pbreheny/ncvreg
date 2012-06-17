plot.cv.ncvreg <- function(x, log.l=TRUE, ...)
  { 
    N <- nrow(x$E)
    e.sd <- apply(x$E,2,sd)
    L <- x$cve - e.sd/sqrt(N)
    U <- x$cve + e.sd/sqrt(N)
    l <- x$lambda
    if (log.l)
      {
        l <- log(l)
        xlab <- expression(log(lambda))
      }
    else xlab <- expression(lambda)

    plot.args = list(x=l, y=x$cve, ylim=range(c(L,U)), xlab=xlab, ylab="Cross-validation error", type="n", xlim=rev(range(l)))
    new.args = list(...)
    if (length(new.args)) plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    abline(v=l[x$min],lty=2,lwd=.5)
    arrows(x0=l,x1=l,y0=L,y1=U,code=3,angle=90,col="gray80",length=.05)
    points(l,x$cve,col="red",pch=19,cex=.5)
  }
