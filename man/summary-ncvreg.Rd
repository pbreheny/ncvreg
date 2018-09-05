\name{summary.ncvreg}
\alias{summary.ncvreg}
\alias{summary.ncvsurv}
\alias{print.summary.ncvreg}
\title{Summary method for ncvreg objects}
\description{Inferential summaries for \code{ncvreg} and \code{ncvsurv}
  objects based on local marginal false discovery rates.}
\usage{
\method{summary}{ncvreg}(object, lambda, which, number, cutoff, ...) 
\method{print}{summary.ncvreg}(x, digits, ...)
}
\arguments{
  \item{object}{An \code{ncvreg} or \code{ncvsurv} object.}
  \item{lambda}{The regularization parameter value at which inference
    should be reported.}
  \item{which}{Alternatively, \code{lambda} may be specified by index;
    \code{which=10} means: report inference for the 10th value of
    \code{lambda} along the regularization path.  If both \code{lambda}
    and \code{which} are specified, \code{lambda} takes precedence.}
  \item{number}{By default, \code{summary} will provide an inferential
    summary for each variable that has been selected (i.e.,  each
    variable with a nonzero coefficient).  Specifying \code{number=5},
    for example, means that the summary table will include the 5
    features with the lowest mfdr values, regardless of whether they
    were selected.  To see all features, \code{number=Inf}.}
  \item{cutoff}{Alternatively, specifying for example \code{cutoff=0.3}
    will report inference for all features with mfdr under 30\%.  If
    both \code{number} and \code{cutoff} are specified, the intersection
    between both sets of features is reported.}
  \item{x}{A \code{summary.ncvreg} object.}
  \item{digits}{Number of digits past the decimal point to print out.
    Can be a vector specifying different display digits for each of the
    five non-integer printed values.}
  \item{\dots}{Further arguments; in particular, if you have set
    \code{returnX=FALSE}, you will need to supply \code{X} and \code{y}
    in order to calculate local mFDRs.}
}
\value{
  \code{summary.ncvreg} and \code{summary.ncvsurv} produce object with
  S3 class \code{summary.ncvreg}.  The class has its own print method
  and contains the following list elements:
  \item{penalty}{The penalty used by \code{ncvreg} or \code{ncvsurv}.}
  \item{model}{Either \code{"linear"}, \code{"logistic"}, or
    \code{"Cox"}.}
  \item{n}{Number of instances.}
  \item{p}{Number of regression coefficients (not including the
    intercept).}
  \item{lambda}{The \code{lambda} value at which inference is being
    reported.}
  \item{nvars}{The number of nonzero coefficients (again, not including
    the intercept) at that value of \code{lambda}.}
  \item{table}{A table containing estimates, normalized test statistics
    (z), and an estimate of the local mfdr for each coefficient.  The
    mfdr may be loosely interpreted, in an empirical Bayes sense, as the
    probability that the given feature is null.}
  \item{unpen.table}{If there are any unpenalized coefficients, a
    separate inferential summary is given for them.  Currently, this is
    based on \code{lm}/\code{glm}/\code{coxph} using the penalized
    coefficients to provide an offset.  This is useful and more or less
    accurate, but not ideal; we hope to improve the inferential methods
    for unpenalized variables in the future.}}
%\references{COMING SOON}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{ncvreg}}, \code{\link{cv.ncvreg}},
  \code{\link{plot.cv.ncvreg}}} 
\examples{
# Linear regression --------------------------------------------------
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y)
summary(fit, lambda=0.08)

# Logistic regression ------------------------------------------------
data(Heart)
fit <- ncvreg(Heart$X, Heart$y, family="binomial")
summary(fit, lambda=0.05)

# Cox regression -----------------------------------------------------
data(Lung)
fit <- ncvsurv(Lung$X, Lung$y)
summary(fit, lambda=0.1)

# Options ------------------------------------------------------------
fit <- ncvreg(Heart$X, Heart$y, family="binomial")
summary(fit, lambda=0.08, number=3)
summary(fit, lambda=0.08, number=Inf)
summary(fit, lambda=0.08, cutoff=0.5)


fit <- ncvreg(Heart$X, Heart$y, family="binomial", returnX=FALSE)
summary(fit, X=Heart$X, y=Heart$y, lambda=0.08)
}
