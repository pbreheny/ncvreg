#' Plot coefficients from a ncvreg object
#' 
#' Produces a plot of the coefficient paths for a fitted \code{ncvreg} object.
#' 
#' @param x Fitted \code{"ncvreg"} model.
#' @param alpha Controls alpha-blending, helpful when the number of covariates
#' is large.  Default is alpha=1.
#' @param log.l Should horizontal axis be on the log scale?  Default is FALSE.
#' @param shade Should nonconvex region be shaded?  Default is TRUE.
#' @param col Vector of colors for coefficient lines.  By default, evenly
#' spaced colors are selected automatically.
#' @param \dots Other graphical parameters to \code{plot}
#' @author Patrick Breheny
#' @seealso \code{\link{ncvreg}}
#' @references Breheny P and Huang J. (2011) Coordinate descentalgorithms for
#' nonconvex penalized regression, with applications to biological feature
#' selection.  \emph{Annals of Applied Statistics}, \strong{5}: 232-253.
#' c("\\Sexpr[results=rd]{tools:::Rd_expr_doi(\"#1\")}",
#' "10.1214/10-AOAS388")\Sexpr{tools:::Rd_expr_doi("10.1214/10-AOAS388")}
#' @examples
#' 
#' data(Prostate)
#' 
#' fit <- ncvreg(Prostate$X, Prostate$y)
#' plot(fit)
#' plot(fit, col="black")
#' plot(fit, log=TRUE)
#' @export

plot.ncvreg <- function(x, alpha=1, log.l=FALSE, shade=TRUE, col, ...) {
  if (length(x$lambda) == 1) stop("Object was fit with only a single lambda value; there is no path to plot", call.=FALSE)
  YY <- if (length(x$penalty.factor)==nrow(x$beta)) coef(x) else coef(x)[-1, , drop=FALSE]
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
    polygon(x=c(l1,l2,l2,l1), y=c(plot.args$ylim[1], plot.args$ylim[1], plot.args$ylim[2], plot.args$ylim[2]), col="gray85", border=FALSE)
  }
  
  if (missing(col)) {
    col <- grDevices::hcl(h=seq(15, 375, len=max(4, p+1)), l=60, c=150, alpha=alpha)
    col <- if (p==2) col[c(1,3)] else col[1:p]
  } else {
    col <- col[ind]
  }
  line.args <- list(col=col, lwd=1+2*exp(-p/20), lty=1)
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- l
  line.args$y <- t(Y)
  do.call("matlines", line.args)
  
  abline(h=0)
}
