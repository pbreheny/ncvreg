\name{perm.ncvreg}
\alias{perm.ncvreg}
\title{Permutation fitting for ncvreg}
\description{Fits multiple penalized regression models in which the
  outcome is randomly permuted, thereby allowing estimation of the false
  inclusion rate.}
\usage{
perm.ncvreg(X, y, ..., N=10, seed, trace=FALSE)
}
\arguments{
  \item{X}{The design matrix, without an intercept, as in
    \code{ncvreg}.}
  \item{y}{The response vector, as in \code{ncvreg}.}
  \item{...}{Additional arguments to \code{ncvreg}.}
  \item{N}{The number of permutation replications.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{trace}{If set to TRUE, perm.ncvreg will inform the user of its
    progress by announcing the beginning of each permutation fit.
    Default is FALSE.}
  }
\details{
  The function fits a penalized regression model to the actual data,
  then repeats the process \code{N} times with a permuted version of the
  response vector.  This allows estimation of the expected number of
  variables included by chance for each value of \code{lambda}.  The
  ratio of this expected quantity to the number of selected variables
  using the actual (non-permuted) response is called the false inclusion
  rate (FIR).}
\value{
  An object with S3 class \code{"perm.ncvreg"} containing:
  \item{EF}{The number of variables selected at each value of
    \code{lambda}, averaged over the permutation fits.}
  \item{S}{The actual number of selected variables for the non-permuted
    data.}
  \item{FIR}{The estimated false inclusion rate (\code{EF/S}).}
  \item{fit}{The fitted \code{ncvreg} object for the original
    (non-permuted) data.} 
  \item{loss}{The loss/deviance for each value of \code{lambda},
    averaged over the permutation fits.  This is an estimate of the
    explanatory power of the model under null conditions, and can be
    used to adjust the loss of the fitted model in a manner akin to
    the idea of an adjusted R-squared in classical regression.}
  }
%\references{None yet}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{ncvreg}}, \code{\link{plot.fir}}, \code{\link{fir}}}
\examples{
## Linear regression
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa
pmfit <- perm.ncvreg(X, y)

par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)
plot(pmfit$fit)

## Logistic regression
data(heart)
X <- as.matrix(heart[,1:9])
y <- heart$chd
pmfit <- perm.ncvreg(X, y, family="binomial")

par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)
plot(pmfit$fit)
}
