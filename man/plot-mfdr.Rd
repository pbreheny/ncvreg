\name{plot.mfdr}
\alias{plot.mfdr}
\title{Plot marginal false discovery rate curves}
\description{Plot marginal false discovery rate curves from an \code{"mfdr"} or
  \code{"perm.ncvreg"} object.}
\usage{
\method{plot}{mfdr}(x, type=c("mFDR", "EF"), log.l=FALSE, selected=TRUE,
legend=TRUE, ...)
}
\arguments{
  \item{x}{A \code{"perm.ncvreg"} or \code{"mfdr"} object.}
  \item{type}{What to plot on the vertical axis.  \code{mFDR} plots the
    marginal false discovery rate; \code{EF} plots the expected number of false
    discoveries along with the actual number of variables included in the
    model.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    FALSE.}
  \item{selected}{If \code{TRUE} (the default), places an axis on top of
    the plot denoting the number of variables in the model (i.e., that
    have a nonzero regression coefficient) at that value of
    \code{lambda}.}
  \item{legend}{For \code{type="EF"} plots, draw a legend to indicate
    which line is for the actual selections and which line is for the
    expected number of false discoveries?  Default is \code{TRUE}.}
  \item{\dots}{Other graphical parameters to pass to \code{plot}}
}
%\details{}
\references{
    Breheny P (to appear). Marginal false discovery rates for penalized
    regression models. Biostatistics.
}
\author{Patrick Breheny}
\seealso{\code{\link{mfdr}}, \code{\link{perm.ncvreg}}}
\examples{
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y)

obj <- mfdr(fit)
obj[1:10,]

# Some plotting options
plot(obj)
plot(obj, type="EF")
plot(obj, log=TRUE)


# Comparison with perm.ncvreg
op <- par(mfrow=c(2,2))
plot(obj)
plot(obj, type="EF")
pmfit <- perm.ncvreg(Prostate$X, Prostate$y)
plot(pmfit)
plot(pmfit, type="EF")
par(op)
}
