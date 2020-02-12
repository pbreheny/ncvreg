#' Standardizes a design matrix
#' 
#' The function \code{std} accepts a design matrix and returns a standardized version of that matrix (i.e., each column will have mean 0 and mean sum of squares equal to 1).
#' 
#'   This function centers and scales each column of \code{X} so that
#'   \deqn{\sum_{i=1}^n x_{ij}=0}{sum(X[,j])=0}
#'   and
#'   \deqn{n^{-1} \sum_{i=1}^n x_{ij}^2 = 1}{mean(X[,j]^2)=1}
#'   for all j.  This is usually not necessary to call directly, as \code{ncvreg} internally standardizes the design matrix, but inspection of the standardized design matrix can sometimes be useful.  This differs from the base R function \code{\link[base]{scale}} in two ways:
#'   
#'   1. `scale` uses the sample standard deviation \code{sqrt(sum(x^2)/(n-1))}, while \code{std} uses the root-mean-square (population) standard deviation \code{sqrt(mean(sum(x^2)))}
#'   2. `std` is faster.
#' 
#' @param X   A matrix (or object that can be coerced to a matrix, such as a data frame or numeric vector).
#' 
#' @return The standardized design matrix, with the following attribues:
#'   * `center`, `scale`: mean and standard deviation used to scale the columns
#'   * `nonsingular`: A vector indicating which columns of the original design matrix were able to be standardized (constant columns cannot be standardized to have a standard deviation of 1)
#' 
#' @examples 
#' X <- matrix(rnorm(50), 10, 5)
#' S <- std(X)
#' apply(S, 2, sum)
#' apply(S, 2, function(x) mean(x^2))

std <- function(X) {
  if (typeof(X) == 'integer') storage.mode(X) <- 'double'
  if (!inherits(X, "matrix")) {
    if (is.numeric(X)) {
      X <- matrix(as.double(X), ncol=1)
    } else {
      tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
  }
  STD <- .Call("standardize", X)
  dimnames(STD[[1]]) <- dimnames(X)
  ns <- which(STD[[3]] > 1e-6)
  if (length(ns) == ncol(X)) {
    val <- STD[[1]]
  } else {
    val <- STD[[1]][, ns, drop=FALSE]
  }
  attr(val, "center") <- STD[[2]]
  attr(val, "scale") <- STD[[3]]
  attr(val, "nonsingular") <- ns
  val
}
