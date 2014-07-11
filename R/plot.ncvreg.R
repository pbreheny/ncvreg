plot.ncvreg <- function(x, alpha=1, log.l=FALSE, shade=TRUE, ...) {
  YY <- if (length(x$penalty.factor)==nrow(x$beta)) coef(x) else coef(x)[-1,,drop=FALSE]
  penalized <- which(x$penalty.factor!=0)
  nonzero <- which(apply(abs(YY), 1, sum)!=0)
  ind <- intersect(penalized, nonzero)
  Y <- YY[ind, , drop=FALSE]
  p <- nrow(Y)
  l <- x$lambda
  
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)
  plot.args <- list(x=l, y=1:length(l), ylim=range(Y), xlab=xlab, ylab="", type="n", xlim=rev(range(l)), las=1)
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  if (!is.element("ylab", names(new.args))) mtext(expression(hat(beta)), side=2, cex=par("cex"), line=3, las=1)
  
  if (shade & !is.null(x$convex.min)) {
    l1 <- l[x$convex.min]
    l2 <- min(l)
    polygon(x=c(l1,l2,l2,l1),y=c(plot.args$ylim[1],plot.args$ylim[1],plot.args$ylim[2],plot.args$ylim[2]),col="gray85",border=FALSE) 
  }
  
  cols <- hcl(h=seq(15, 375, len=max(4, p+1)), l=60, c=150, alpha=alpha)
  cols <- if (p==2) cols[c(1,3)] else cols[1:p]  
  line.args <- list(col=cols, lwd=1+2*exp(-p/20), lty=1)
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- l
  line.args$y <- t(Y)
  do.call("matlines",line.args)
  
  abline(h=0)
}
