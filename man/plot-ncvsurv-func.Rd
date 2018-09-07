\name{plot.ncvsurv.func}
\alias{plot.ncvsurv.func}
\title{Plot survival curve for ncvsurv model}
\description{Plot survival curve for a model that has been fit using
  \code{ncvsurv} followed by a prediction of the survival function using
  \code{predict.ncvsurv}}
\usage{
\method{plot}{ncvsurv.func}(x, alpha=1, ...)
}
\arguments{
  \item{x}{A \code{'ncvsurv.func'} object, which is returned by
    \code{predict.ncvsurv} if \code{type='survival'} is specified.  See
    examples.}
  \item{alpha}{Controls alpha-blending (i.e., transparency).  Useful if many
    overlapping lines are present.}
  \item{\dots}{Other graphical parameters to pass to \code{plot}}
}
%\details{}
%\references{None yet}
\author{Patrick Breheny}
\seealso{\code{\link{ncvsurv}}, \code{\link{predict.ncvsurv}}}
\examples{
data(Lung)
X <- Lung$X
y <- Lung$y

fit <- ncvsurv(X, y)

# A single survival curve
S <- predict(fit, X[1,], type='survival', lambda=.15)
plot(S, xlim=c(0,200))

# Lots of survival curves
S <- predict(fit, X, type='survival', lambda=.08)
plot(S, xlim=c(0,200), alpha=0.3)
}
