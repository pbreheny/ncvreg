plot.ncvsurv.func <- function(x, alpha=1, ...) {
  time <- attr(x, 'time')
  if (length(x) > 1) {
    Y <- sapply(x, function(f) f(time))
    n <- ncol(Y)
  } else {
    Y <- x(time)
    n <- 1
  }

  plot.args <- list(x=1, y=1, xlim=range(time), ylim=range(Y), xlab='Time', ylab='Pr(Survival)', type="n", las=1)
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)

  cols <- hcl(h=seq(15, 375, len=max(4, n+1)), l=60, c=150, alpha=alpha)
  cols <- if (n==2) cols[c(1,3)] else cols[1:n]
  line.args <- list(col=cols, lwd=1+2*exp(-n/20), lty=1, type='s')
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- time
  line.args$y <- Y
  do.call("matlines",line.args)
}
