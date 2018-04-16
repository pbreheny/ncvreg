\name{predict.ncvsurv}
\alias{predict.ncvsurv}
\alias{coef.ncvsurv}
\title{Model predictions based on a fitted "ncvsurv" object.}
\description{
Similar to other predict methods, this function returns predictions from
a fitted \code{"ncvsurv"} object.
}
\usage{
\method{predict}{ncvsurv}(object, X, type=c("link", "response", "survival",
"median", "coefficients", "vars", "nvars"), lambda,
which=1:length(object$lambda), ...)
}
\arguments{
  \item{object}{Fitted \code{"ncvsurv"} model object.}
  \item{X}{Matrix of values at which predictions are to be made.  Not
    used for \code{type="coefficients"} or for some of the \code{type}
    settings in \code{predict}.}
  \item{lambda}{Values of the regularization parameter \code{lambda} at
    which predictions are requested.  For values of \code{lambda} not in
    the sequence of fitted models, linear interpolation is used.}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    predictions are required.  By default, all indices are returned.  If
    \code{lambda} is specified, this will override \code{which}.}
  \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the risk (i.e., exp(link));
    \code{"survival"} returns the estimated survival function; \code{"median"}
    estimates median survival times.  The other options are all identical to
    their \code{ncvreg} counterparts: \code{"coefficients"} returns the
    coefficients; \code{"vars"} returns a list containing the indices and names
    of the nonzero variables at each value of \code{lambda}; \code{"nvars"}
    returns the number of nonzero coefficients at each value of \code{lambda}.}
  \item{\dots}{Not used.}}
\details{Estimation of baseline survival function conditional on the estimated
  values of \code{beta} is carried out according to the method described in
  Chapter 4.3 of Kalbfleish and Prentice.  In particular, it agrees exactly
  the results returned by \code{survfit.coxph(..., type='kalbfleisch-prentice')}
  in the \code{survival} package.}
\value{The object returned depends on type.}
\references{
  \itemize{
    \item Breheny P and Huang J (2011). Coordinate descent algorithms for
      nonconvex penalized regression, with applications to biological feature
      selection. Ann. Appl. Statist., 5: 232-253.
    \item Kalbfleish JD and Prentice RL (2002). The Statistical Analysis of
      Failure Time Data, 2nd edition. Wiley.
  }
}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{ncvsurv}}}
\examples{
data(Lung)
X <- Lung$X
y <- Lung$y

fit <- ncvsurv(X,y)
coef(fit, lambda=0.05)
head(predict(fit, X, type="link", lambda=0.05))
head(predict(fit, X, type="response", lambda=0.05))

# Survival function
S <- predict(fit, X[1,], type="survival", lambda=0.05)
S(100)
S <- predict(fit, X, type="survival", lambda=0.05)
plot(S, xlim=c(0,200))

# Medians
predict(fit, X[1,], type="median", lambda=0.05)
M <- predict(fit, X, type="median")
M[1:10, 1:10]

# Nonzero coefficients
predict(fit, type="vars", lambda=c(0.1, 0.01))
predict(fit, type="nvars", lambda=c(0.1, 0.01))
}
