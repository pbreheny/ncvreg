\name{AUC.cv.ncvsurv}
\alias{AUC}
\alias{AUC.cv.ncvsurv}
\title{Calculates AUC for cv.ncvsurv objects}
\description{Calculates the cross-validated AUC (concordance) from a
  "cv.ncvsurv" object.}
\usage{
\method{AUC}{cv.ncvsurv}(obj, ...)
}
\arguments{
\item{obj}{A \code{cv.ncvsurv} object.  You must run \code{cv.ncvsurv}  with the option \code{returnY=TRUE} in order for \code{AUC}  to work.}

\item{\dots}{For S3 method compatibility?}
}
\details{
  The area under the curve (AUC), or equivalently, the concordance
  statistic (C), is calculated according to the procedure outlined in
  the reference below.  This calls the \code{survConcordance} function
  in the \code{survival} package, except the cross-validated linear
  predictors are used to guard against overfitting.  Thus, the values
  returned by \code{AUC.cv.ncvsurv} will be lower than those
  you would obtain with \code{survConcordance} if you fit the full
  (unpenalized) model.}
\references{van Houwelingen H, Putter H (2011). Dynamic Prediction in
  Clinical Survival Analysis.  CRC Press.}
\author{Patrick Breheny, Brandon Butcher, and Lawrence Hunsicker}
\seealso{\code{\link{cv.ncvsurv}}, \code{survConcordance}}
\examples{
data(Lung)
X <- Lung$X
y <- Lung$y

cvfit <- cv.ncvsurv(X, y, returnY=TRUE)
head(AUC(cvfit))
lam <- cvfit$lambda
plot(lam, AUC(cvfit), xlim=rev(range(lam)), lwd=3, type='l',
     las=1, xlab=expression(lambda), ylab='AUC')
}
