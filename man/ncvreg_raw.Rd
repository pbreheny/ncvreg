\name{ncvreg_raw}
\alias{ncvreg_raw}
\title{Fit an MCP- or SCAD-penalized regression path on raw input data}
\description{Fit coefficients paths for MCP- or SCAD-penalized
  regression models on raw input data over a grid of values for the regularization
  parameter lambda.  Currently supports linear regression, with
  option for an additional L2 penalty.}
\usage{
ncvreg_raw(X, y, family=c("gaussian", "binomial", "poisson"),
penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, 3),
alpha=1, lambda.min=ifelse(n>p,.0001,.01), nlambda=100, lambda, eps=1e-4,
max.iter=10000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(X)),
warn=TRUE, returnX=FALSE, intercept=TRUE, ...)
}
\arguments{
  \item{X}{The design matrix, without an intercept.  \code{ncvreg_raw}
  uses X as is without standardization.}
  \item{y}{The response vector.}
  \item{family}{Either "gaussian", "binomial", or "poisson", depending
    on the response. Currently only "gaussian" is supported.}
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
    lambda.max.  Default is .0001 if the number of observations is larger
    than the number of covariates and .01 otherwise.}
  \item{nlambda}{The number of lambda values.  Default is 100.}
  \item{lambda}{A user-specified sequence of lambda values.  By default,
    a sequence of values of length \code{nlambda} is computed, equally
    spaced on the log scale.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    RMSD for the change in linear predictors for any coefficient is less
    than \code{eps}.  Default is \code{1e-4}.}
  \item{max.iter}{Maximum number of iterations (total across entire path).
    Default is 1000.}
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
    always in the model without shrinkage. Intercept is not penalized 
    when one is included by setting \code{intercept} to TRUE.}
  \item{warn}{Return warning messages for failures to converge and model
    saturation?  Default is TRUE.}
  \item{returnX}{Return the standardized design matrix?  Default is
    FALSE.}
  \item{intercept}{Include an intercept or not? Default is TRUE.}
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
  for \code{"gaussian"}.\cr
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
    (\code{"gaussian"}) or negative log-likelihood (\code{"binomial"}
    and \code{"poisson"}) of the fitted model at each value of
    \code{lambda}.}
  \item{penalty.factor}{Same as above.}
}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Congrui Yi <congrui-yi@uiowa.edu>}
\seealso{\code{\link{ncvreg}}, \code{\link{plot.ncvreg}}, \code{\link{cv.ncvreg}}}
\examples{
## Linear regression
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

par(mfrow=c(2,1))
fit <- ncvreg_raw(X,y)
plot(fit,main="With Intercept")
fit <- ncvreg_raw(X,y,intercept = FALSE)
plot(fit,main="Without Intercept")
}
