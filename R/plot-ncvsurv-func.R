#' Plot survival curve for ncvsurv model
#' 
#' Plot survival curve for a model that has been fit using \code{ncvsurv}
#' followed by a prediction of the survival function using
#' \code{predict.ncvsurv}
#' 
#' 
#' @param x A \code{'ncvsurv.func'} object, which is returned by
#' \code{predict.ncvsurv} if \code{type='survival'} is specified.  See
#' examples.
#' @param alpha Controls alpha-blending (i.e., transparency).  Useful if many
#' overlapping lines are present.
#' @param \dots Other graphical parameters to pass to \code{plot}
#' @author Patrick Breheny
#' @seealso \code{\link{ncvsurv}}, \code{\link{predict.ncvsurv}}
#' @examples
#' 
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#' 
#' fit <- ncvsurv(X, y)
#' 
#' # A single survival curve
#' S <- predict(fit, X[1,], type='survival', lambda=.15)
#' plot(S, xlim=c(0,200))
#' 
#' # Lots of survival curves
#' S <- predict(fit, X, type='survival', lambda=.08)
#' plot(S, xlim=c(0,200), alpha=0.3)
#' @export

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

  cols <- grDevices::hcl(h=seq(15, 375, len=max(4, n+1)), l=60, c=150, alpha=alpha)
  cols <- if (n==2) cols[c(1,3)] else cols[1:n]
  line.args <- list(col=cols, lwd=1+2*exp(-n/20), lty=1, type='s')
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- time
  line.args$y <- Y
  do.call("matlines", line.args)
}
