#' Standardizes a design matrix
#' 
#' Accepts a design matrix and returns a standardized version of that matrix
#' (i.e., each column will have mean 0 and mean sum of squares equal to 1).
#' 
#' This function centers and scales each column of `X` so that
#' \deqn{\sum_{i=1}^n x_{ij}=0}
#' and
#' \deqn{n^{-1} \sum_{i=1}^n x_{ij}^2 = 1}
#' for all j. This is usually not necessary to call directly, as **ncvreg** internally
#' standardizes the design matrix, but inspection of the standardized design matrix
#' can sometimes be useful. This differs from the base R function [scale()] in two ways:
#'   
#'   1. `scale()` uses the sample standard deviation `sqrt(sum(x^2)/(n-1))`, while `std()` uses the root-mean-square standard deviation `sqrt(mean(sum(x^2))` without the \eqn{n/(n-1)} correction
#'   2. `std` is faster.
#' 
#' @param X     A matrix (or object that can be coerced to a matrix, such as a data frame or numeric vector).
#' @param Xnew  Optional. If supplied, `X` must be the output of `std()` and `Xnew` is to be standardized in the same way. See examples for why this might be useful.
#' 
#' @returns The standardized design matrix, with the following attribues:
#'   \item{center, scale}{mean and standard deviation used to scale the columns}
#'   \item{nonsingular}{A vector indicating which columns of the original design matrix were able to be standardized (constant columns cannot be standardized to have a standard deviation of 1)}
#' 
#' @examples 
#' data(Prostate)
#' S <- std(Prostate$X)
#' apply(S, 2, sum)
#' apply(S, 2, function(x) mean(x^2))
#' 
#' # Standardizing new observations
#' X1 <- Prostate$X[1:90,]
#' X2 <- Prostate$X[91:97,]
#' S <- std(X1)
#' head(std(S, X2))
#' # Useful if you fit to a standardized X, but then get new obs:
#' y <- Prostate$y[1:90]
#' fit <- ncvreg(S, y)
#' predict(fit, std(S, X2), lambda=0.1)
#' # Same as
#' predict(ncvreg(X1, y), X2, lambda=0.1)
#' @export

std <- function(X, Xnew) {
  if (missing(Xnew)) {
    if (typeof(X) == 'integer') storage.mode(X) <- 'double'
    if (!inherits(X, "matrix")) {
      if (is.numeric(X)) {
        X <- matrix(as.double(X), ncol=1)
      } else {
        tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
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
  } else {
    attrx <- names(attributes(X))
    if (!('center' %in% attrx & 'scale' %in% attrx & 'nonsingular' %in% attrx)) {
      stop('If supplying both X and Xnew, must run std(X) first and pass the result to std', call.=FALSE)
    }
    val <- scale(Xnew[,attr(X, 'nonsingular')], attr(X, 'center'), attr(X, 'scale'))
    attributes(val)[c('nonsingular', 'center', 'scale')] <- attributes(X)[c('nonsingular', 'center', 'scale')]
    attributes(val)[c('scaled:center', 'scaled:scale')] <- NULL
  }
  val
}
