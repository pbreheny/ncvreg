#' Pipe NCVREG Results into Confidence Intervals and Significance Testing
#'
#' This function processes the output of a cross-validated `ncvreg` fit, extracting the estimated coefficients, computing standard errors, and producing confidence intervals along with significance testing.
#'
#' @param cv_fit An object of class `cv.ncvreg` containing the cross-validated fit. Must include components `fit`, `lambda.min`, and other relevant information.
#' @param alpha A numeric value specifying the significance level for confidence intervals. Default is 0.05.
#' @param X (Optional) The design matrix. If NULL, it is extracted from `cv_fit`.
#' @param y (Optional) The response vector. If NULL, it is extracted from `cv_fit`.
#'
#' @return A data.frame containing the following columns:
#' \describe{
#'   \item{estimate}{The estimated coefficients (divided by rescale).}
#'   \item{lower}{The lower bounds of the confidence intervals.}
#'   \item{upper}{The upper bounds of the confidence intervals.}
#'   \item{significance}{Adjusted p-values for significance testing, corrected using the Benjamini-Hochberg method.}
#'   \item{original_pvals}{The original p-values before adjustment.}
#' }
#'
#' @details The function assumes the input `cv_fit` contains all necessary components, including the design matrix `X` and the response vector `y`. The function recalculates the response mean, residuals, and applies significance testing using the Benjamini-Hochberg correction.
#'
#' @seealso \code{\link[ncvreg]{ncvreg}}, \code{\link[ncvreg]{cv.ncvreg}}
#'
#' @examples
#' \dontrun{
#' library(ncvreg)
#' X <- matrix(rnorm(100*20), 100, 20)
#' y <- rnorm(100)
#' cvfit <- cv.ncvreg(X, y, family="gaussian")
#' result <- pipe_ncvreg(cvfit)
#' print(result)
#' }
#'
#' @export
pipe_ncvreg <- function(cv_fit, alpha = 0.05, X = NULL, y = NULL) {
  
  if (class(cv_fit) != "cv.ncvreg") {
    stop("cv_fit must be an opbject of class cv.ncvreg.")
  }
 
  # Check if cv_fit$fit contains X and y, or if they are supplied
  if ((is.null(cv_fit$fit$X) & is.null(X)) & ((is.null(cv_fit$fit$y) & is.null(y)))) {
    stop("Please supply X and y, or specify returnX = TRUE and returnY = TRUE in your call to cv.ncvreg.")
  } else if (((is.null(cv_fit$fit$X) & is.null(X)) | ((is.null(cv_fit$fit$y) & is.null(y))))) {
    if (((is.null(cv_fit$fit$X) & is.null(X))) {
      stop("Please supply X, or specify returnX = TRUE in your call to cv.ncvreg.")
    } else if ((is.null(cv_fit$fit$y) & is.null(y))) {
      stop("Please supply Y, or specify returnY = TRUE in your call to cv.ncvreg.")
    }
  } else {
    if ((!is.null(X) & !is.null(cv_fit$fit$X)) & (!is.null(y) & !is.null(cv_fit$fit$y))) {
      message("Both cv_fit$fit and user-supplied X and y are provided; using cv_fit$fit's X and y.")
    } else if (!is.null(X) & !is.null(cv_fit$fit$X)) {
      # Message indicating preference for X and y from cv_fit$fit
      message("Both cv_fit$fit and user-supplied X are provided; using cv_fit$fit's X.")
    } else if (!is.null(y) & !is.null(cv_fit$fit$y)) {
      # Message indicating preference for X and y from cv_fit$fit
      message("Both cv_fit$fit and user-supplied y are provided; using cv_fit$fit's y.")
    }
  }
  
  n <- cv_fit$fit$n
  X <- cv_fit$fit$X
  y <- cv_fit$fit$y
  
  bh_lambda <- coef(cv_fit$fit, lambda = cv_fit$lambda.min)
  rescale <- attr(X, "scale")
  bh_lambda <- bh_lambda[-1] * rescale
  intercept <- mean(y - as.numeric(X %*% bh_lambda))
  p <- length(bh_lambda)
  yhat <- intercept + as.numeric(X %*% bh_lambda)
  
  
  s <- bh_lambda != 0
  sh_lh <- sum(s)
  sigma_h <- sqrt((n - sh_lh)^(-1) * sum((y - yhat)^2))
  
  partial_residuals <- (y - yhat) + (X * matrix(bh_lambda, nrow = n, ncol = p, byrow=TRUE))
  b_bar <- (1/n)*colSums(X * partial_residuals)
  
  if (sh_lh > 0) {
    Xs <- X[,s,drop =FALSE]
    q_sh_default <- diag(n) - Xs %*% solve(t(Xs) %*% Xs, tol = 1e-12) %*% t(Xs)
  } else {
    q_sh_default <- diag(n)
  }
  
  ses <- numeric(p)
  for (j in 1:p) {
    
    if (!s[j]) {
      s_j <- s
      sh_j <- sh_lh
      q_sh <- q_sh_default
    } else {
      s_j <- s
      s_j[j] <- FALSE
      sh_j <- sum(s_j)
      
      if (sh_j > 0) {
        Xsj <- X[,s_j,drop=FALSE]
        q_sh <- diag(n) - Xsj %*% solve(t(Xsj) %*% Xsj, tol = 1e-12) %*% t(Xsj)
      } else {
        q_sh <- diag(n)
      }
      
    }
    
    ses[j] <- sigma_h * sqrt((t(X[,j,drop=FALSE]) %*% q_sh %*% X[,j,drop=FALSE])^(-1))
    
  }
  
  ts <- b_bar / ses
  ps <- 2 * (1 - pnorm(abs(ts)))
  qs <- p.adjust(ps, method = "BH")
  widths <- abs(qnorm(alpha / 2)) * ses
  
  data.frame(estimate = b_bar / rescale, lower = (b_bar - widths) / rescale, upper = (b_bar + widths) / rescale, significance = qs, original_pvals = ps)
  
}