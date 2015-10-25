\name{auc.cv.ncvsurv}
\alias{auc}
\alias{auc.cv.ncvsurv}
\title{Calculates AUC for cv.ncvsurv objects}
\description{Calculates the cross-validated AUC (concordance) from a
  "cv.ncvsurv" object.}
\usage{
\method{auc}{cv.ncvsurv}(obj, ...)
}
\arguments{
  \item{obj}{A \code{cv.ncvsurv} object.  You must run \code{cv.ncvsurv}
  with the option \code{returnY=TRUE} in order for \code{auc} to work.}
  \item{\dots}{For S3 method compatibility.}
}
\details{
  The area under the curve (AUC), or equivalently, the concordance
  statistic (C), is calculated according to the procedure outlined in
  the reference below.  This is the same procedure used by the
  \code{survConcordance} function in the \code{survival} package, except
  for the fact that the cross-validated linear predictors are used to
  guard against overfitting.  Thus, the values returned by
  \code{auc.cv.ncvsurv} will be somewhat lower than those you would
  obtain with \code{survConcordance} if you fit the full (unpenalized)
  model.
  }
\references{van Houwelingen H, Putter H (2011). Dynamic Prediction in
  Clinical Survival Analysis.  CRC Press.}
\author{Brandon Butcher; imported into the ncvreg package by Patrick
  Breheny}
\seealso{\code{\link{cv.ncvsurv}}}
\examples{
data(Lung)
X <- Lung$X
y <- Lung$y

cvfit <- cv.ncvsurv(X, y, returnY=TRUE)
head(auc(cvfit))

ll <- log(cvfit$fit$lambda)
plot(ll, auc(cvfit), xlim=rev(range(ll)), lwd=3, type='l',
     xlab=expression(log(lambda)), ylab='AUC')
}
