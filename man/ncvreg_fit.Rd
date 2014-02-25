\name{ncvreg_fit}
\alias{ncvreg_fit}
\title{Internal C function for fitting ncvreg models}
\description{This function is a wrapper for the calls to C carried out
  by \code{ncvreg}.  This function performs no checking, processing or
  standardization, so use \code{ncvreg} instead, unless you really know
  what you're doing.}
\usage{
ncvreg_fit(X, y, family=c("gaussian", "binomial", "poisson"),
penalty=c("MCP", "SCAD", "lasso"), gamma=3, alpha=1, lambda, eps=.001,
max.iter=1000, dfmax=p+1, penalty.factor=rep(1, ncol(X)), warn=TRUE)
}
\arguments{
  \item{X}{The design matrix.  Unlike in \code{\link{ncvreg}},
    \code{ncvreg_fit} does not standardize the data.  For Gaussian
    responses, no intercept is included.  At least for now, an intercept
    is still included automatically for binomial and Poisson responses,
    as it is not clear (to me) that removing the intercept here would
    make sense.} 
  \item{y}{The response vector.}
  \item{family}{Either "gaussian", "binomial", or "poisson", depending
    on the response.}
  \item{penalty}{The penalty to be applied to the model.  Either "MCP"
    (the default), "SCAD", or "lasso".}
  \item{gamma}{The tuning parameter of the MCP/SCAD penalty, as in
    \code{\link{ncvreg}}.}
  \item{alpha}{Tuning parameter which controls the relative
    contributions from the MCP/SCAD penalty and the ridge, 
    or L2 penalty, as in \code{\link{ncvreg}}.}
  \item{lambda}{A user-specified sequence of lambda values.  Unlike in
    \code{\link{ncvreg}}, this must be specified directly in
    \code{ncvreg_fit}.}
  \item{eps}{Convergence threshhold, as in \code{\link{ncvreg}}.}
  \item{max.iter}{Maximum number of iterations.  Default is 1000.}
  \item{dfmax}{Upper bound for the number of nonzero coefficients, as in
    \code{\link{ncvreg}}.}
  \item{penalty.factor}{A multiplicative factor for the penalty applied
    to each coefficient, as in \code{\link{ncvreg}}.}
  \item{warn}{Return warning messages for failures to converge and model
    saturation?  Default is TRUE.}
}
\details{
  \code{ncvreg_fit} is supplied as a separate function in case
  developers wish to embed ncvreg's internal algorithms in a larger
  procedure.  It should not be called directly unless you know exactly
  what you are doing.  In particular, no standardization or processing
  of \code{X} and \code{y} are carried out, and the output will not work
  with any of other functions in the package such as
  \code{\link{plot.ncvreg}} or \code{\link{predict.ncvreg}}.}
\value{
  A list with components:
  \item{beta}{The fitted matrix of coefficients.  The number of rows is
    equal to the number of coefficients, and the number of columns is
    equal to the length of \code{lambda}.}
  \item{loss}{A vector containing either the residual sum of squares
    (\code{"gaussian"}) or negative log-likelihood (\code{"binomial"} or
    \code{"poisson"}) of the fitted model at each value of
    \code{lambda}.}
  \item{iter}{A vector of length \code{nlambda} containing the number
    of iterations until convergence at each value of \code{lambda}.}
  \item{lambda}{The sequence of regularization parameter values in the
    path.}
  }
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{ncvreg}}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

## These two results are NOT the same;
## No standardization is being done in the latter
ncvreg(X, y, lambda=c(0.5, 0.1, 0.05))$beta
ncvreg_fit(X, y, lambda=c(0.5, 0.1, 0.05))$beta
}
