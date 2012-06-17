\name{predict.ncvreg}
\alias{coef.ncvreg}
\alias{predict.ncvreg}
\title{Model predictions based on a fitted "ncvreg" object.}
\usage{
  \method{predict}{ncvreg} (object, X, lambda,
    which = 1:length(object$lambda),
    type = c("link", "response", "class", "coefficients"),
    ...)
}
\arguments{
  \item{object}{Fitted \code{"ncvreg"} model object.}

  \item{X}{Matrix of values at which predictions are to be
  made.  Not used for \code{type="coefficients"}}

  \item{lambda}{Value of tuning parameter at which
  predictions are requested.  By default, returns all
  values obtained during fit.  Overrides \code{which}.}

  \item{which}{Indices of the penalty parameter
  \code{lambda} at which predictions are required.  By
  default, all indices are returned.}

  \item{type}{Type of prediction: \code{"link"} returns the
  linear predictors; \code{"response"} gives the fitted
  values; \code{"class"} returns the binomial outcome with
  the highest probability; \code{"coefficients"} returns
  the coefficients.}

  \item{...}{Not used.}
}
\value{
  The object returned depends on type.
}
\description{
  Similar to other predict methods, this function returns
  predictions from a fitted \code{"ncvreg"} object.
}
\examples{
data(heart)
X <- as.matrix(heart[,1:9])
y <- heart$chd

par(mfrow=c(2,2))
fit <- ncvreg(X,y,family="binomial")
coef(fit)
predict(fit,X,type="link")
predict(fit,X,type="response")
predict(fit,X,type="class")
}
\author{
  Patrick Breheny <patrick.breheny@uky.edu>
}
\references{
  Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with
  applications to biological feature selection. Ann. Appl.
  Statist., 5: 232-253.
}
\seealso{
  \code{ncvreg}
}
\keyword{models}
\keyword{regression}

