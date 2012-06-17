\name{ncvreg-package}
\alias{ncvreg-package}
\docType{package}
\title{Regularization paths for SCAD- and MCP-penalized regression models}
\description{Efficient algorithms for fitting regularization paths for
  linear or logistic regression models penalized by MCP or SCAD, with
  optional additional L2 penalty ("Mnet").}
\details{
  \tabular{ll}{
    Package: \tab ncvreg\cr
    Type: \tab Package\cr
    Version: \tab 2.3\cr
    Date: \tab 2011-05-06\cr
    License: \tab GPL-2\cr}
  Accepts a design matrix \code{X} and vector of responses \code{y},
  produces the regularization path over a grid of values for the tuning
  parameter \code{lambda}.  Also provides methods for plotting and for
  determining locally convex regions of the coefficients paths.}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\keyword{models}
\keyword{regression}
\keyword{package}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

fit <- ncvreg(X,y)
plot(fit)
}
