plot.ncvreg <- function(x, alpha=1, log.l=FALSE, shade=TRUE, ...)
  {
    zeros <- which(apply(abs(x$beta),1,sum)==0)
    beta <- x$beta[-c(1,zeros),,drop=FALSE]
    p <- nrow(beta)
    l <- x$lambda
    
    if (log.l)
      {
        l <- log(l)
        xlab <- expression(log(lambda))
      }
    else xlab <- expression(lambda)
    plot.args <- list(x=l, y=1:length(l), ylim=range(beta), xlab=xlab, ylab=expression(hat(beta)), type="n", xlim=rev(range(l)))
    new.args <- list(...)
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)

    if (shade & !is.null(x$convex.min))
      {
        l1 <- l[x$convex.min]
        l2 <- min(l)
        polygon(x=c(l1,l2,l2,l1),y=c(plot.args$ylim[1],plot.args$ylim[1],plot.args$ylim[2],plot.args$ylim[2]),col="gray85",border=FALSE) 
      }
    
    line.args <- list(col=hcl(h=seq(0,360,len=(p+1)),l=70,c=100,alpha=alpha)[1:p],lwd=1+1.2^(-p/20),lty=1)
    if (length(new.args)) line.args[names(new.args)] <- new.args
    line.args$x <- l
    line.args$y <- t(beta)
    do.call("matlines",line.args)
    
    abline(h=0)
  }
