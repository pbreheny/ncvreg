\name{plot.ncvreg}
\alias{plot.ncvreg}
\title{Plot coefficients from a "ncvreg" object}
\description{Produces a plot of the coefficient paths for a fitted
  \code{"ncvreg"} object.}
\usage{
\method{plot}{ncvreg}(x, alpha=1, log.l=FALSE, shade=TRUE, ...)
}
\arguments{
  \item{x}{Fitted \code{"ncvreg"} model.}
  \item{alpha}{Controls alpha-blending, helpful when the number of
    covariates is large.  Default is alpha=1.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    FALSE.} 
  \item{shade}{Should nonconvex region be shaded?  Default is TRUE.}
  \item{\dots}{Other graphical parameters to \code{plot}}
}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{ncvreg}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

fit <- ncvreg(X,y)
plot(fit)
plot(fit,col="black")
plot(fit,log=TRUE)
}
\keyword{models}
\keyword{regression}
