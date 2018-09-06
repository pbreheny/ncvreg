\name{cv.ncvreg}
\alias{cv.ncvreg}
\title{Cross-validation for ncvreg}
\description{Performs k-fold cross validation for MCP- or SCAD-penalized
  regression models over a grid of values for the regularization
  parameter lambda.}
\usage{
cv.ncvreg(X, y, ..., cluster, nfolds=10, seed, fold, returnY=FALSE,
trace=FALSE) }
\arguments{
  \item{X}{The design matrix, without an intercept, as in
    \code{ncvreg}.}
  \item{y}{The response vector, as in \code{ncvreg}.}
  \item{...}{Additional arguments to \code{ncvreg}.}
  \item{cluster}{\code{cv.ncvreg} can be run in parallel across a
    cluster using the \code{parallel} package.  The cluster must be set
    up in advance using the \code{makeCluster} function from that
    pacakge.  The cluster must then be passed to \code{cv.ncvreg} (see
    example).}
  \item{nfolds}{The number of cross-validation folds.  Default is 10.}
  \item{fold}{Which fold each observation belongs to.  By default the
    observations are randomly assigned by \code{cv.ncvreg}.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{returnY}{Should \code{cv.ncvreg} return the fitted values from
    the cross-validation folds?  Default is FALSE; if TRUE, this will
    return a matrix in which the element for row i, column j is the
    fitted value for observation i from the fold in which observation i
    was excluded from the fit, at the jth value of lambda.}
  \item{trace}{If set to TRUE, cv.ncvreg will inform the user of its
    progress by announcing the beginning of each CV fold.  Default is
    FALSE.}
  }
\details{
  The function calls \code{ncvreg} \code{nfolds} times, each time
  leaving out 1/\code{nfolds} of the data.  The cross-validation
  error is based on the residual sum of squares when
  \code{family="gaussian"} and the deviance when
  \code{family="binomial"} or \code{family="poisson"}.
  \href{http://pbreheny.github.io/ncvreg/articles/web/models.html}{More
  details are available on the ncvreg website}.
  
  For \code{family="binomial"} models, the cross-validation fold
  assignments are balanced across the 0/1 outcomes, so that each fold
  has the same proportion of 0/1 outcomes (or as close to the same
  proportion as it is possible to achieve if cases do not divide
  evenly).
}
\value{
  An object with S3 class \code{"cv.ncvreg"} containing:
  \describe{
    \item{cve}{The error for each value of \code{lambda}, averaged
      across the cross-validation folds.}
    \item{cvse}{The estimated standard error associated with each value of
      for \code{cve}.}
    \item{lambda}{The sequence of regularization parameter values along
      which the cross-validation error was calculated.}
    \item{fit}{The fitted \code{ncvreg} object for the whole data.}
    \item{min}{The index of \code{lambda} corresponding to
      \code{lambda.min}.}
    \item{lambda.min}{The value of \code{lambda} with the minimum
      cross-validation error.}
    \item{null.dev}{The deviance for the intercept-only model.}
    \item{Bias}{The estimated bias of the minimum cross-validation error,
      as in Tibshirani RJ and Tibshirani R (2009),
      "A Bias Correction for the Minimum Error Rate in Cross-Validation",
      Ann. Appl. Stat. 3:822-829.}
    \item{pe}{If \code{family="binomial"}, the cross-validation prediction
      error for each value of \code{lambda}.}
    \item{Y}{If \code{returnY=TRUE}, the matrix of cross-validated fitted
      values (see above).}
  }
}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny; Grant Brown helped with the parallelization
  support}
\seealso{\code{\link{ncvreg}}, \code{\link{plot.cv.ncvreg}},
  \code{\link{summary.cv.ncvreg}}}
\examples{
data(Prostate)

cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
plot(cvfit)
summary(cvfit)

fit <- cvfit$fit
plot(fit)
beta <- fit$beta[,cvfit$min]

## requires loading the parallel package
\dontrun{
library(parallel)
X <- Prostate$X
y <- Prostate$y
cl <- makeCluster(4)
cvfit <- cv.ncvreg(X, y, cluster=cl, nfolds=length(y))}
}
