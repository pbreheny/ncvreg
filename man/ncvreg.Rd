\name{ncvreg}
\alias{ncvreg}
\title{Fit an MCP- or SCAD-penalized regression path}
\description{Fit coefficients paths for MCP- or SCAD-penalized
  regression models over a grid of values for the regularization
  parameter lambda.  Fits linear and logistic regression models, with
  option for an additional L2 penalty.}
\usage{
ncvreg(X, y, family=c("gaussian","binomial"), penalty=c("MCP", "SCAD",
"lasso"), gamma=switch(penalty, SCAD=3.7, 3), alpha=1,
lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=.001,
max.iter=1000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(X)),
warn=TRUE, returnX=FALSE, ...)
}
\arguments{
  \item{X}{The design matrix, without an intercept.  \code{ncvreg}
  standardizes the data and includes an intercept by default.}
  \item{y}{The response vector.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{penalty}{The penalty to be applied to the model.  Either "MCP"
    (the default), "SCAD", or "lasso".}
  \item{gamma}{The tuning parameter of the MCP/SCAD penalty (see
    details).  Default is 3 for MCP and 3.7 for SCAD.}
  \item{alpha}{Tuning parameter for the Mnet estimator which controls
    the relative contributions from the MCP/SCAD penalty and the ridge,
    or L2 penalty.  \code{alpha=1} is equivalent to MCP/SCAD penalty,
    while \code{alpha=0} would be equivalent to ridge regression.
    However, \code{alpha=0} is not supported; \code{alpha} may be
    arbitrarily small, but not exactly 0.}
  \item{lambda.min}{The smallest value for lambda, as a fraction of
    lambda.max.  Default is .001 if the number of observations is larger
    than the number of covariates and .05 otherwise.}
  \item{nlambda}{The number of lambda values.  Default is 100.}
  \item{lambda}{A user-specified sequence of lambda values.  By default,
    a sequence of values of length \code{nlambda} is computed, equally
    spaced on the log scale.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{.001}.}
  \item{max.iter}{Maximum number of iterations.  Default is 1000.x}
  \item{convex}{Calculate index for which objective function ceases to
    be locally convex?  Default is TRUE.}
  \item{dfmax}{Upper bound for the number of nonzero coefficients.
    Default is no upper bound.  However, for large data sets,
    computational burden may be heavy for models with a large number of
    nonzero coefficients.}
  \item{penalty.factor}{A multiplicative factor for the penalty applied
    to each coefficient.  If supplied, \code{penalty.factor} must be a
    numeric vector of length equal to the number of columns of
    \code{X}.  The purpose of \code{penalty.factor} is to apply
    differential penalization if some coefficients are thought to be
    more likely than others to be in the model.  In particular,
    \code{penalty.factor} can be 0, in which case the coefficient is
    always in the model without shrinkage.}
  \item{warn}{Return warning messages for failures to converge and model
    saturation?  Default is TRUE.}
  \item{returnX}{Return the standardized design matrix?  Default is
    FALSE.}
  \item{...}{Not used.}
}
\details{
  The sequence of models indexed by the regularization parameter
  \code{lambda} is fit using a coordinate descent algorithm.  For
  logistic regression models, some care is taken to avoid model
  saturation; the algorithm may exit early in this setting.  The
  objective function is defined to be
  \deqn{\frac{1}{2n}\textrm{RSS} + \textrm{penalty}}{RSS/(2*n) +
    penalty} 
  for \code{"gaussian"} and
  \deqn{-\frac{1}{n}\ell + \textrm{penalty}}{-(1/n) loglik + penalty}
  for \code{"binomial"}, where the likelihood is from a traditional
  generalized linear model for the log-odds of an event.\cr
  \cr
  This algorithm is stable, very efficient, and generally converges
  quite rapidly to the solution.  For logistic regression, adaptive
  rescaling (see reference) is used.\cr
  \cr
  The convexity diagnostics rely on a fine covering of
  (lambda.min,lambda.max); choosing a low value of \code{nlambda} may
  produce unreliable results.}
\value{
  An object with S3 class \code{"ncvreg"} containing:
  \item{beta}{The fitted matrix of coefficients.  The number of rows is
    equal to the number of coefficients, and the number of columns is
    equal to \code{nlambda}.}
  \item{iter}{A vector of length \code{nlambda} containing the number
    of iterations until convergence at each value of \code{lambda}.}
  \item{lambda}{The sequence of regularization parameter values in the
    path.}
  \item{penalty}{Same as above.}
  \item{family}{Same as above.}
  \item{gamma}{Same as above.}
  \item{alpha}{Same as above.}
  \item{convex.min}{The last index for which the objective function is
    locally convex.  The smallest value of lambda for which the
    objective function is convex is therefore \code{lambda[convex.min]},
    with corresponding coefficients \code{beta[,convex.min]}.}
  \item{loss}{A vector containing either the residual sum of squares
    (\code{"gaussian"}) or negative log-likelihood (\code{"binomial"})
    of the fitted model at each value of \code{lambda}.}}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny <patrick-breheny@uiow.edu>}
\seealso{\code{plot.ncvreg}, \code{cv.ncvreg}}
\examples{
## Linear regression
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

par(mfrow=c(2,2))
fit <- ncvreg(X,y)
plot(fit,main=expression(paste(gamma,"=",3)))
fit <- ncvreg(X,y,gamma=10)
plot(fit,main=expression(paste(gamma,"=",10)))
fit <- ncvreg(X,y,gamma=1.5)
plot(fit,main=expression(paste(gamma,"=",1.5)))
fit <- ncvreg(X,y,penalty="SCAD")
plot(fit,main=expression(paste("SCAD, ",gamma,"=",3)))

par(mfrow=c(2,2))
fit <- ncvreg(X,y)
plot(fit,main=expression(paste(alpha,"=",1)))
fit <- ncvreg(X,y,alpha=0.9)
plot(fit,main=expression(paste(alpha,"=",0.9)))
fit <- ncvreg(X,y,alpha=0.5)
plot(fit,main=expression(paste(alpha,"=",0.5)))
fit <- ncvreg(X,y,alpha=0.1)
plot(fit,main=expression(paste(alpha,"=",0.1)))

par(mfrow=c(2,2))
fit <- ncvreg(X,y)
plot(fir(fit))            ## Independence approximation
plot(fir(fit), type="EF") ## Independence approximation
perm.fit <- perm.ncvreg(X,y)
plot(perm.fit)
plot(perm.fit, type="EF")

## Logistic regression
data(heart)
X <- as.matrix(heart[,1:9])
y <- heart$chd

par(mfrow=c(2,2))
fit <- ncvreg(X,y,family="binomial")
plot(fit,main=expression(paste(gamma,"=",3)))
fit <- ncvreg(X,y,family="binomial",gamma=10)
plot(fit,main=expression(paste(gamma,"=",10)))
fit <- ncvreg(X,y,family="binomial",gamma=1.5)
plot(fit,main=expression(paste(gamma,"=",1.5)))
fit <- ncvreg(X,y,family="binomial",penalty="SCAD")
plot(fit,main=expression(paste("SCAD, ",gamma,"=",3)))

par(mfrow=c(2,2))
fit <- ncvreg(X,y,family="binomial")
plot(fit,main=expression(paste(alpha,"=",1)))
fit <- ncvreg(X,y,family="binomial",alpha=0.9)
plot(fit,main=expression(paste(alpha,"=",0.9)))
fit <- ncvreg(X,y,family="binomial",alpha=0.5)
plot(fit,main=expression(paste(alpha,"=",0.5)))
fit <- ncvreg(X,y,family="binomial",alpha=0.1)
plot(fit,main=expression(paste(alpha,"=",0.1)))
}
\keyword{models}
\keyword{regression}
