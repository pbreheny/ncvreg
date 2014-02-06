\name{plot.fir}
\alias{plot.fir}
\title{Plot false inclusion rate curves}
\description{Plot false inclusion rate curves from an \code{"fir"} or
  \code{"perm.ncvreg"} object.}
\usage{
\method{plot}{fir}(x, type=c("FIR", "EF"), log.l=FALSE, legend=TRUE,
...)
}
\arguments{
  \item{x}{A \code{"perm.ncvreg"} or \code{"fir"} object.}
  \item{type}{What to plot on the vertical axis.  \code{FIR} plots the
    false inclusion rate; \code{EF} plots the expected number of false
    inclusions along with the actual number of variables included in the
    model.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    FALSE.}
  \item{legend}{For \code{type="EF"} plots, draw a legend to indicate
    which line is for the actual selections and which line is for the
    expected number of false inclusions?  Default is \code{TRUE}.}
  \item{\dots}{Other graphical parameters to pass to \code{plot}}
}
%\details{}
%\references{None yet}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{fir}}, \code{\link{perm.ncvreg}}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa
fit <- ncvreg(X, y)

f <- fir(fit)
cbind(EF=f$EF, S=f$S, FIR=f$FIR)[1:20,]

## Comparison with perm.ncvreg
par(mfrow=c(2,2))
plot(f)
plot(f, type="EF")
pmfit <- perm.ncvreg(X, y)
plot(pmfit)
plot(pmfit, type="EF")
## Note that fir() is more conservative
}
