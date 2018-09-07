\name{mfdr}
\alias{mfdr}
\title{Marginal false discovery rates}
\description{Estimates the marginal false discovery rate (mFDR) of a penalized
  regression model.}
\usage{
mfdr(fit, X)
}
\arguments{
  \item{fit}{An \code{ncvreg} or \code{ncvsurv} object.}
  \item{X}{The model matrix corresponding to \code{fit}.  This is not
    necessary for linear regression, but in logistic and Cox regression,
    the mFDR depends on X.  It is not necessary to supply \code{X} if it
    is already contained in \code{fit}; i.e., if
    \code{ncvreg}/\code{ncvsurv} was run with \code{returnX=TRUE}.}
  }
\details{
  The function estimates the marginal false discovery rate (mFDR) for a
  penalized regression model.  The estimate tends to be accurate in most
  settings, but will be slightly conservative if predictors are highly
  correlated.  For an alternative way of estimating the mFDR, typically more
  accurate in highly correlated cases,  see \code{\link{perm.ncvreg}}.}
\value{
  An object with S3 class \code{mfdr} inheriting from \code{data.frame} and
  containing:
  \item{EF}{The number of variables selected at each value of
    \code{lambda}, averaged over the permutation fits.}
  \item{S}{The actual number of selected variables for the non-permuted
    data.}
  \item{mFDR}{The estimated marginal false discovery rate (\code{EF/S}).}
  }
%\references{None yet}
\author{Patrick Breheny and Ryan Miller}
\seealso{\code{\link{ncvreg}}, \code{\link{ncvsurv}}, \code{\link{plot.mfdr}},
  \code{\link{perm.ncvreg}}}
\examples{

# Linear regression --------------------------------
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y)

obj <- mfdr(fit)
obj[1:10,]

# Comparison with perm.ncvreg
op <- par(mfrow=c(2,2))
plot(obj)
plot(obj, type="EF")
pmfit <- perm.ncvreg(Prostate$X, Prostate$y)
plot(pmfit)
plot(pmfit, type="EF")
par(op)

# Logistic regression ------------------------------
data(Heart)
fit <- ncvreg(Heart$X, Heart$y, family="binomial")
obj <- mfdr(fit)
head(obj)
op <- par(mfrow=c(1,2))
plot(obj)
plot(obj, type="EF")
par(op)

# Cox regression -----------------------------------
data(Lung)
fit <- ncvsurv(Lung$X, Lung$y)
obj <- mfdr(fit)
head(obj)
op <- par(mfrow=c(1,2))
plot(obj)
plot(obj, type="EF")
par(op)
}
