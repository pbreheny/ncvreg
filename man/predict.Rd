\name{predict.ncvreg}
\alias{predict.ncvreg}
\alias{coef.ncvreg}
\title{Model predictions based on a fitted ncvreg object.}
\description{
Similar to other predict methods, this function returns predictions from
a fitted \code{ncvreg} object.
}
\usage{
\method{predict}{ncvreg}(object, X, type=c("link", "response", "class",
"coefficients", "vars", "nvars"), lambda, which=1:length(object$lambda),
...)
\method{coef}{ncvreg}(object, lambda, which=1:length(object$lambda),
drop=TRUE, ...)
}
\arguments{
  \item{object}{Fitted \code{ncvreg} model object.}
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
    linear predictors; \code{"response"} gives the fitted values;
    \code{"class"} returns the binomial outcome with the highest
    probability; \code{"coefficients"} returns the coefficients;
    \code{"vars"} returns a list containing the indices and names of the
    nonzero variables at each value of \code{lambda}; \code{"nvars"}
    returns the number of nonzero coefficients at each value of
    \code{lambda}.}
  \item{drop}{If coefficients for a single value of \code{lambda} are to
    be returned, reduce dimensions to a vector?  Setting
    \code{drop=FALSE} returns a 1-column matrix.}
  \item{\dots}{Not used.}}
\value{The object returned depends on type.}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny}
\seealso{\code{\link{ncvreg}}}
\examples{
data(Heart)

fit <- ncvreg(Heart$X, Heart$y, family="binomial")
coef(fit, lambda=0.05)
head(predict(fit, Heart$X, type="link", lambda=0.05))
head(predict(fit, Heart$X, type="response", lambda=0.05))
head(predict(fit, Heart$X, type="class", lambda=0.05))
predict(fit, type="vars", lambda=c(0.05, 0.01))
predict(fit, type="nvars", lambda=c(0.05, 0.01))
}
