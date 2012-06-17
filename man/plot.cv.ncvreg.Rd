\name{plot.cv.ncvreg}
\alias{plot.cv.ncvreg}
\title{Plots the cross-validation curve from a "cv.ncvreg" object}
\description{Plots the cross-validation curve from a "cv.ncvreg"
  object, along with standard error bars.}
\usage{
\method{plot}{cv.ncvreg}(x, log.l=TRUE, ...)
}
\arguments{
  \item{x}{A \code{"cv.ncvreg"} object.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    TRUE.} 
  \item{\dots}{Other graphical parameters to \code{plot}}
}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{ncvreg}, \code{cv.ncvreg}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

cvfit <- cv.ncvreg(X,y)
plot(cvfit)
}
\keyword{models}
\keyword{regression}
