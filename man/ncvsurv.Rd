\name{ncvsurv}
\alias{ncvsurv}
\title{Fit an MCP- or SCAD-penalized survival model}
\description{Fit coefficients paths for MCP- or SCAD-penalized
  Cox regression models over a grid of values for the regularization
  parameter lambda, with option for an additional L2 penalty.}
\usage{
ncvsurv(X, y, penalty=c("MCP", "SCAD", "lasso"),
gamma=switch(penalty, SCAD=3.7, 3), alpha=1,
lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=1e-4,
max.iter=10000, convex=TRUE, dfmax=p, penalty.factor=rep(1, ncol(X)),
warn=TRUE, returnX, ...)
}
\arguments{
  \item{X}{The design matrix of predictor values.  \code{ncvsurv}
    standardizes the data prior to fitting.}
  \item{y}{The time-to-event outcome, as a two-column matrix or
    \code{\link[survival]{Surv}} object.  The first column should be time
      on study (follow up time); the second column should be a binary
      variable with 1 indicating that the event has occurred and 0
      indicating (right) censoring.}
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
    RMSD for the change in linear predictors for any coefficient is less
    than \code{eps}.  Default is \code{1e-4}.}
  \item{max.iter}{Maximum number of iterations (total across entire
    path). Default is 1000.}
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
    always in the model without any penalization/shrinkage.}
  \item{warn}{Return warning messages for failures to converge and model
    saturation?  Default is TRUE.}
  \item{returnX}{Return the standardized design matrix along with the
    fit?  By default, this option is turned on if X is under 100 MB, but
    turned off for larger matrices to preserve memory.  Note that
    certain methods, such as \code{\link{summary.ncvsurv}} require
    access to the design matrix and may not be able to run if
    \code{returnX=FALSE}.}
  \item{...}{Not used.}
}
\details{
  The sequence of models indexed by the regularization parameter
  \code{lambda} is fit using a coordinate descent algorithm.  In order
  to accomplish this, the second derivative (Hessian) of the Cox partial
  log-likelihood is diagonalized (see references for details).  The
  objective function is defined to be
  \deqn{Q(\beta|X, y) = \frac{1}{n} L(\beta|X, y) +
    P_\lambda(\beta)}{Q(\beta|X, y) = (1/n)*L(\beta|X, y) +
    P(\beta, \lambda),}
  where the loss function L is the deviance (-2 times the partial
  log-likelihood) from the Cox regression mode.
  \href{http://pbreheny.github.io/ncvreg/articles/web/models.html}{See
    here for more details}.

  Presently, ties are not handled by \code{ncvsurv} in a particularly
  sophisticated manner.  This will be improved upon in a future release
  of \code{ncvreg}.
}
\value{
  An object with S3 class \code{"ncvsurv"} containing:
  \describe{
    \item{beta}{The fitted matrix of coefficients.  The number of rows
      is equal to the number of coefficients, and the number of columns
      is equal to \code{nlambda}.}
    \item{iter}{A vector of length \code{nlambda} containing the number
      of iterations until convergence at each value of \code{lambda}.}
    \item{lambda}{The sequence of regularization parameter values in the
      path.}
    \item{penalty}{Same as above.}
    \item{model}{Same as above.}
    \item{gamma}{Same as above.}
    \item{alpha}{Same as above.}
    \item{convex.min}{The last index for which the objective function is
      locally convex.  The smallest value of lambda for which the
      objective function is convex is therefore
      \code{lambda[convex.min]}, with corresponding coefficients
      \code{beta[,convex.min]}.}
    \item{loss}{The deviance of the fitted model at each value of
      \code{lambda}.}
    \item{penalty.factor}{Same as above.}
    \item{n}{The number of observations.}
  }

  For Cox models, the following objects are also returned (and are
  necessary to estimate baseline survival conditonal on the estimated
  regression coefficients), all of which are ordered by time on study.
  I.e., the ith row of \code{W} does not correspond to the ith row of
  \code{X}):

  \describe{
    \item{W}{Matrix of \code{exp(beta)} values for each subject over all
      \code{lambda} values.}
    \item{time}{Times on study.}
    \item{fail}{Failure event indicator.}
  }

  Additionally, if \code{returnX=TRUE}, the object will also contain

  \describe{
    \item{X}{The standardized design matrix.}
  }
}
\references{
  \itemize{
    \item Breheny P and Huang J. (2011) Coordinate descentalgorithms
    for nonconvex penalized regression, with applications to biological
    feature selection.  \emph{Annals of Applied Statistics}, \strong{5}:
    232-253.
    \url{myweb.uiowa.edu/pbreheny/publications/Breheny2011.pdf}

    \item Simon N, Friedman JH, Hastie T, and Tibshirani R. (2011)
    Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent.  \emph{Journal of Statistical Software},
    \strong{39}: 1-13.
    \url{http://www.jstatsoft.org/v39/i05}
  }
}
\author{Patrick Breheny}
\seealso{\code{\link{plot.ncvreg}}, \code{\link{cv.ncvsurv}}}
\examples{
data(Lung)
X <- Lung$X
y <- Lung$y

op <- par(mfrow=c(2,2))
fit <- ncvsurv(X, y)
plot(fit, main=expression(paste(gamma,"=",3)))
fit <- ncvsurv(X, y, gamma=10)
plot(fit, main=expression(paste(gamma,"=",10)))
fit <- ncvsurv(X, y, gamma=1.5)
plot(fit, main=expression(paste(gamma,"=",1.5)))
fit <- ncvsurv(X, y, penalty="SCAD")
plot(fit, main=expression(paste("SCAD, ",gamma,"=",3)))
par(op)

fit <- ncvsurv(X,y)
ll <- log(fit$lambda)
op <- par(mfrow=c(2,1))
plot(ll, BIC(fit), type="l", xlim=rev(range(ll)))
lam <- fit$lambda[which.min(BIC(fit))]
b <- coef(fit, lambda=lam)
b[b!=0]
plot(fit)
abline(v=lam)
par(op)

S <- predict(fit, X, type='survival', lambda=lam)
plot(S, xlim=c(0,200))
}
