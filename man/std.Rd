\name{std}
\alias{std}
\title{Standardizes a design matrix}
\description{The function \code{std} accepts a design matrix and returns
  a standardized version of that matrix (i.e., each column will have
  mean 0 and mean sum of squares equal to 1).}
\usage{
std(X)
}
\arguments{
  \item{X}{A matrix (or object that can be coerced to a matrix,
    such as a data frame).}
}
\details{
  This function centers and scales each column of \code{X} so that
  \deqn{\sum_{i=1}^n x_{ij}=0}{sum(X[,j])=0}
  and
  \deqn{n^{-1} \sum_{i=1}^n x_{ij}^2 = 1}{mean(X[,j]^2)=1}
  for all j.  This is usually not necessary to call directly, as \code{ncvreg}
  internally standardizes the design matrix, but inspection of the
  standardized design matrix can sometimes be useful.  This differs from
  the base R function \code{\link[base]{scale}} in two ways: (1) \code{scale}
  uses the sample standard deviation \code{sqrt(sum(x^2)/(n-1))}, while
  \code{std} uses the root-mean-square, or population, standard
  deviation \code{sqrt(mean(sum(x^2)))}, and (2) \code{std} is faster.
}
\value{
  The standardized design matrix, with attribues "center" and "scale"
  corresponding to the mean and (population) standard deviation used to
  scale the columns.
}
\author{Patrick Breheny}
\examples{
X <- matrix(rnorm(50), 10, 5)
S <- std(X)
apply(S, 2, sum)
apply(S, 2, function(x) mean(x^2))
}
