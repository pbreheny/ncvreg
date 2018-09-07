\name{permres}
\alias{permres}
\alias{permres.ncvreg}
\title{Permute residuals for a fitted ncvreg model}
\description{Fits multiple penalized regression models in which the
  residuals are randomly permuted, thereby allowing estimation of the
  marginal false discovery rate.}
\usage{
permres(fit, \dots)
\method{permres}{ncvreg}(fit, lambda, N=10, seed, trace=FALSE, \dots)
}
\arguments{
  \item{fit}{A fitted ncvreg model, as produced by
    \code{\link{ncvreg}()}. To use with \code{permres}, the model must
    be fit using the \code{returnX=TRUE} option.}
  \item{lambda}{The regularization parameter to use for estimating
    residuals. Unlike \code{\link{perm.ncvreg}}, \code{permres}
    calculates EF and mFDR for a specific \code{lambda} value, not an
    entire path.  As a result, it runs much faster.}
  \item{N}{The number of permutation replications.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{trace}{If set to TRUE, perm.ncvreg will inform the user of its
    progress by announcing the beginning of each permutation
    fit. Default is FALSE.}
  \item{\dots}{Not used.}
}
\details{
  The function fits a penalized regression model to the actual data,
  then repeats the process \code{N} times with a permuted version of the
  response vector.  This allows estimation of the expected number of
  variables included by chance for each value of \code{lambda}.  The
  ratio of this expected quantity to the number of selected variables
  using the actual (non-permuted) response is called the marginal false
  discovery rate (mFDR).}
\value{
  A list with the following components:
  \item{EF}{The number of variables selected at each value of
    \code{lambda}, averaged over the permutation fits.}
  \item{S}{The actual number of selected variables for the non-permuted
    data.}
  \item{mFDR}{The estimated marginal false discovery rate (\code{EF/S}).}
  \item{loss}{The loss/deviance, averaged over the permutation fits.  This is an
    estimate of the explanatory power of the model under null conditions, and
    can be used to adjust the loss of the fitted model in a manner akin to
    the idea of an adjusted R-squared in classical regression.}
  }
%\references{None yet}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{ncvreg}}, \code{\link{mfdr}}, \code{\link{perm.ncvreg}}}
\examples{
data(Prostate)
fit <- ncvreg(Prostate$X, Prostate$y, N=50)
permres(fit, lambda=0.15)
}
