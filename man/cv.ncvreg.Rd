\name{cv.ncvreg}
\alias{cv.ncvreg}
\title{Cross-validation for ncvreg}
\description{Performs k-fold cross validation for MCP- or SCAD-penalized
  regression models over a grid of values for the regularization
  parameter lambda.}
\usage{
cv.ncvreg(X, y, family=c("gaussian","binomial"), alpha=1,
lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, nfolds=10, seed,
trace=FALSE, ...)
}
\arguments{
  \item{X}{The design matrix, without an intercept.  \code{ncvreg}
  standardizes the data and includes an intercept by default.}
  \item{y}{The response vector.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{alpha}{Tuning parameter for the Mnet estimator which controls
    the relative contributions from the MCP/SCAD penalty and the ridge,
    or L2 penalty.  \code{alpha=1} is equivalent to MCP/SCAD penalty,
    while \code{alpha=0} would be equivalent to ridge regression.
    However, \code{alpha=0} is not supported; \code{alpha} may be
    arbitrarily small, but not exactly 0.}
  \item{lambda.min}{The smallest value for lambda, as a fraction of
    lambda.max.  Default is .001 if the number of observations is larger
    than the number of covariates and .05 otherwise.}
  \item{nlambda}{The number of lambda values.  Default is 100.}
  \item{lambda}{A user-specified sequence of lambda values.  By default,
    a sequence of values of length \code{nlambda} is computed, equally
    spaced on the log scale.}
  \item{nfolds}{The number of cross-validation folds.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{trace}{If set to TRUE, cv.ncvreg will inform the user of its
    progress by announcing the beginning of each CV fold.  Default is
    FALSE.}
  \item{...}{Additional arguments to \code{ncvreg}.}}
\details{
  The function calls \code{ncvreg} \code{nfolds} times, each time
  leaving out 1/\code{nfolds} of the data.  The cross-validation
  error is based on the loss function used to fit the model; see
  \code{ncvreg} for details.}
\value{
  An object with S3 class \code{"cv.ncvreg"}.  In addition to the
  \code{ncvreg} fit, the object contains:
  \item{E}{A matrix with \code{nfolds} rows and \code{nlambda}
    columns, containing the cross-validation error for each fit.}
  \item{cve}{The error for each value of \code{lambda}, averaged
    across the cross-validation folds.}
  \item{lambda}{The sequence of regularization parameter values along
    which the cross-validation error was calculated.}
  \item{min}{The index of \code{lambda} corresponding to the minimum
    cross-validation error.}}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{ncvreg}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

cv.fit <- cv.ncvreg(X,y)
plot(cv.fit)

fit <- ncvreg(X,y)
beta <- fit$beta[,cv.fit$min]
}
\keyword{models}
\keyword{regression}
