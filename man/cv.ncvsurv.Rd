\name{cv.ncvsurv}
\alias{cv.ncvsurv}
\title{Cross-validation for ncvsurv}
\description{Performs k-fold cross validation for MCP- or SCAD-penalized
  survival models over a grid of values for the regularization
  parameter lambda.}
\usage{
cv.ncvsurv(X, y, ..., nfolds=10, seed, trace=FALSE
cv.ncvsurv(X, y, ..., nfolds=10, seed, trace=FALSE, events.only=TRUE)
}
\arguments{
  \item{X}{The design matrix, as in \code{ncvsurv}.}
  \item{y}{The response matrix, as in \code{ncvsurv}.}
  \item{...}{Additional arguments to \code{ncvsurv}.}
  \item{nfolds}{The number of cross-validation folds.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{trace}{If set to TRUE, cv.ncvsurv will inform the user of its
    progress by announcing the beginning of each CV fold.  Default is
    FALSE.}
  \item{events.only}{Evaluate cross-validation loss on event times only?
    The default is TRUE.  If set to FALSE, out-of-sample partial
    likelihood deviance is computed for the censored observations as
    well.}
  }
\details{
  The function calls \code{ncvsurv} \code{nfolds} times, each time
  leaving out 1/\code{nfolds} of the data.  Because of the
  semiparametric nature of Cox regression, cross-validation is not
  clearly defined.  \code{cv.ncvsurv} uses the approach proposed by
  Verweij and van Houwelingen (see references), in which the
  contribution to the deviance from observation i comes from subtracting
  the full-data-set partial likelihood from the leave-one-out partial
  likelihood.}
\value{
  An object with S3 class \code{"cv.ncvsurv"} containing:
  \item{cve}{The error for each value of \code{lambda}, averaged
    across the cross-validation folds.}
  \item{cvse}{The estimated standard error associated with each value of
    for \code{cve}.}
  \item{lambda}{The sequence of regularization parameter values along
    which the cross-validation error was calculated.}
  \item{fit}{The fitted \code{ncvsurv} object for the whole data.}
  \item{min}{The index of \code{lambda} corresponding to
    \code{lambda.min}.}
  \item{lambda.min}{The value of \code{lambda} with the minimum
    cross-validation error.}
  \item{null.dev}{The cross-validated deviance for the first model along
    the grid of \code{lambda} (i.e., the cross-validated deviance for
    max(lambda), unless you have supplied your own \code{lambda}
    sequence, in which case this quantity is probably not meaningful).
    Although the actual null deviance can be calculated, it cannot be
    compared with the cross-validated deviance due to the manner in
    which deviance must be calculated for Cox regression models (see
    details).}
  }
\references{
  \itemize{
    \item Breheny P and Huang J. (2011) Coordinate descentalgorithms
    for nonconvex penalized regression, with applications to biological
    feature selection.  \emph{Annals of Applied Statistics}, \strong{5}:
    232-253.
    \url{myweb.uiowa.edu/pbreheny/publications/Breheny2011.pdf}

    \item Verweij PJ and van Houwelingen HC. (1993) Cross-validation in
    survival analysis.  \emph{Statistics in Medicine}, \strong{12}:
    2305-2314.
  }
}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{ncvsurv}, \code{plot.cv.ncvsurv}, \code{summary.cv.ncvsurv}}
\examples{
data(lungCancer)
X <- lungCancer$X
y <- lungCancer$y

cvfit <- cv.ncvsurv(X, y)
summary(cvfit)
plot(cvfit)
plot(cvfit, type="rsq")
}
