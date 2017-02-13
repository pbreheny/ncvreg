\name{mfdr}
\alias{mfdr}
\title{Marginal false discovery rates}
\description{Estimates the marginal false discovery rate (mFDR) of a penalized
  regression model.}
\usage{
mfdr(fit)
}
\arguments{
  \item{fit}{An \code{ncvreg} or \code{ncvsurv} object.}
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
\author{Patrick Breheny <patrick-breheny@uiowa.edu>, with contributions from
  Ryan Miller <ryan-e-miller@uiowa.edu>}
\seealso{\code{\link{ncvreg}}, \code{\link{ncvsurv}}, \code{\link{plot.mfdr}},
  \code{\link{perm.ncvreg}}}
\examples{

#### Linear regression -----------------------------
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa
fit <- ncvreg(X, y)

obj <- mfdr(fit)
obj[1:10,]

# Comparison with perm.ncvreg
par(mfrow=c(2,2))
plot(obj)
plot(obj, type="EF")
pmfit <- perm.ncvreg(X, y)
plot(pmfit)
plot(pmfit, type="EF")
# Note that mfdr() is more conservative

#### Logistic regression ---------------------------
data(heart)
X <- as.matrix(heart[,1:9])
y <- heart$chd
fit <- ncvreg(X, y, family="binomial")
obj <- mfdr(fit, X)
head(obj)
plot(obj)
plot(obj, type="EF")

#### Cox regression --------------------------------
data(Lung)
X <- Lung$X
y <- Lung$y
fit <- ncvsurv(X, y)
obj <- mfdr(fit, X)
head(obj)
plot(obj)
plot(obj, type="EF")
}
