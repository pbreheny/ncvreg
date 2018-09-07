\name{summary.cv.ncvreg}
\alias{summary.cv.ncvreg}
\alias{print.summary.cv.ncvreg}
\title{Summarizing cross-validation-based inference}
\description{Summary method for \code{cv.ncvreg} objects}
\usage{
\method{summary}{cv.ncvreg}(object, ...)
\method{print}{summary.cv.ncvreg}(x, digits, ...)
}
\arguments{
  \item{object}{A \code{"cv.ncvreg"} or \code{"cv.ncvsurv"} object.}
  \item{x}{A \code{"summary.cv.ncvreg"} object.}
  \item{digits}{Number of digits past the decimal point to print out.
    Can be a vector specifying different display digits for each of the
    five non-integer printed values.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
  \code{summary.cv.ncvreg} produces an object with S3 class
  \code{"summary.cv.ncvreg"}.  The class has its own print method and
  contains the following list elements:
  \describe{
    \item{penalty}{The penalty used by \code{ncvreg}.}
    \item{model}{Either \code{"linear"} or \code{"logistic"}, depending on
      the \code{family} option in \code{ncvreg}.}
    \item{n}{Number of observations}
    \item{p}{Number of regression coefficients (not including the
      intercept).}
    \item{min}{The index of \code{lambda} with the smallest
      cross-validation error.}
    \item{lambda}{The sequence of \code{lambda} values used by
      \code{cv.ncvreg}.}
    \item{cve}{Cross-validation error (deviance).}
    \item{r.squared}{Proportion of variance explained by the model, as
      estimated by cross-validation.  For models outside of linear
      regression, the Cox-Snell approach to defining R-squared is used.}
    \item{snr}{Signal to noise ratio, as estimated by cross-validation.}
    \item{sigma}{For linear regression models, the scale parameter
      estimate.}
    \item{pe}{For logistic regression models, the prediction error
      (misclassification error).}
  }
}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny}
\seealso{\code{\link{ncvreg}}, \code{\link{cv.ncvreg}},
  \code{\link{plot.cv.ncvreg}}} 
\examples{
# Linear regression --------------------------------------------------
data(Prostate)
cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
summary(cvfit)

# Logistic regression ------------------------------------------------
data(Heart)
cvfit <- cv.ncvreg(Heart$X, Heart$y, family="binomial")
summary(cvfit)

# Cox regression -----------------------------------------------------
data(Lung)
cvfit <- cv.ncvsurv(Lung$X, Lung$y)
summary(cvfit)
}
