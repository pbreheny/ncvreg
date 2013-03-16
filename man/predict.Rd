\name{predict.ncvreg}
\alias{predict.ncvreg}
\alias{coef.ncvreg}
\title{Model predictions based on a fitted "ncvreg" object.}
\description{
Similar to other predict methods, this function returns predictions from
a fitted \code{"ncvreg"} object.
}
\usage{
\method{predict}{ncvreg}(object, X, lambda,
which=1:length(object$lambda), type=c("link", "response", "class",
"coefficients"), ...) 
\method{coef}{ncvreg}(object, lambda, which=1:length(object$lambda), ...)
}
\arguments{
  \item{object}{Fitted \code{"ncvreg"} model object.}
  \item{X}{Matrix of values at which predictions are to be made.  Not
    used for \code{type="coefficients"}}
  \item{lambda}{Values of the regularization parameter \code{lambda} at
    which predictions are requested.  For values of \code{lambda} not in
    the sequence of fitted models, linear interpolation is used.}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    predictions are required.  By default, all indices are returned.  If
    \code{lambda} is specified, this will override \code{which}.}
  \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the fitted values;
    \code{"class"} returns the binomial outcome with the highest
    probability; \code{"coefficients"} returns the coefficients.}
  \item{\dots}{Not used.}}
\value{The object returned depends on type.}
\references{Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{ncvreg}}
\examples{
data(heart)
X <- as.matrix(heart[,1:9])
y <- heart$chd

fit <- ncvreg(X,y,family="binomial")
coef(fit, lambda=0.05)
predict(fit,X,type="link", lambda=0.05)
predict(fit,X,type="response", lambda=0.05)
predict(fit,X,type="class", lambda=0.05)
}
\keyword{models}
\keyword{regression}
