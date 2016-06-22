\name{fir}
\alias{fir}
\title{False inclusion rates for ncvreg (independence approximation)}
\description{Estimates false inclusion rates (FIR) for penalized
  regression models based on an approximation of independence between
  the predictors.}
\usage{
fir(fit)
}
\arguments{
  \item{fit}{An \code{ncvreg} object.}
  }
\details{
  The function estimates the false inclusion rate (FIR) for a penalized
  regression model.  The calculation is based on an approximation of
  independence between the predictors, and is reasonably accurate in
  near-independent settings.  However, the estimate is conservative when
  predictors are correlated.  For a more accurate estimate of the false
  inclusion rate in the presence of correlated predictors, see
  \code{\link{perm.ncvreg}}.}
\value{
  An object with S3 class \code{fir} inheriting from \code{data.frame} and
  containing:
  \item{EF}{The number of variables selected at each value of
    \code{lambda}, averaged over the permutation fits.}
  \item{S}{The actual number of selected variables for the non-permuted
    data.}
  \item{FIR}{The estimated false inclusion rate (\code{EF/S}).}
  }
%\references{None yet}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{ncvreg}}, \code{\link{plot.fir}}, \code{\link{perm.ncvreg}}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa
fit <- ncvreg(X, y)

f <- fir(fit)
f[1:10,]

## Comparison with perm.ncvreg
par(mfrow=c(2,2))
plot(f)
plot(f, type="EF")
pmfit <- perm.ncvreg(X, y)
plot(pmfit)
plot(pmfit, type="EF")
## Note that fir() is more conservative
}
