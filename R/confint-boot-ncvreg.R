#' Calculating percentile confidence intervals from Hybrid bootstrap draws
#' 
#' Calculates percentile based confidence intervals from Hybrid bootstrap draws
#' obtained from \code{boot_ncvreg}. This produces intervals that are generally
#' conservative on average, especially when n is small or the underlying covarites are
#' especially sparse. However, intervals may undercover when correlation amoung
#' the predictors is high. Although overall coverage is generally achieved,
#' the intervals produced do no attempt to debias and as a result
#' under cover large covariates and over cover covariates near zero.
#'
#' @param object An object of type \code{boot_ncvreg}
#' @param parm A specification of which parameters are to be given confidence 
#' intervals, either a vector of numbers or a vector of names. 
#' If missing, all parameters are considered.
#' @param level The desired significance level to obtain confidence intervals
#' with overall coverage greater than or equal to \code{level}
#' @param ... Not used
#'
#' @return An object with S3 class \code{data.frame} with the columns: \describe{
#' \item{variable}{If the original covariates \code{X} had column names,
#' these are included as the first column.}
#' \item{estimate}{A length \code{ncol(X)} vector of the estimates from the
#' lasso model fit on the original data corresponding to \code{lambda}.}
#' \item{lower}{A length \code{ncol(X)} vector of the lower bounds from the
#' Hybrid bootstrap draws with significance level \code{level}}
#' \item{upper}{A length \code{ncol(X)} vector of the upper bounds from the
#' Hybrid bootstrap draws with significance level \code{level}}}
#' 
#' @examples
#' data(Prostate)
#' boot <- boot_ncvreg(Prostate$X, Prostate$y)
#' confint(boot, level = 0.9)
#' @export
confint.boot_ncvreg <- function(object, parm, level = 0.8, ...) {
  
  if ((missing(object) || class(object) != "boot_ncvreg")) {
    stop("object must be supplied and be of class boot_ncvreg.")
  }
  
  if (!missing(parm)) {
    if (!(is.numeric(parm) || is.character(parm))) {
      stop("Error: 'parm' must be a vector of numbers or a vector of names (characters).")
    }
  }
  
  if (!is.double(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("Error: 'level' must be a double between 0 and 1.")
  }
  
  if (!missing(parm)) {
    draws <- object[["draws"]][,parm,drop=FALSE]
    estimates <- object[["estimates"]][parm]
    variable_names <- names(estimates)
  } else {
    draws <- object[["draws"]] 
    estimates <- object[["estimates"]]
    variable_names <- names(estimates)
  }
  
  any_nas <- any(as.logical(apply(draws, 2, function(x) sum(is.na(x)) > 0)))
  if (any_nas) {
    warning("NAs in draws, this occurs when a varaible is nonsingular (constant) for a given bootstrap draw.")
  }
  
  alpha <- 1 - level
  lowers <- apply(draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  uppers <- apply(draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
    
  ci <- data.frame(
      variable = variable_names,
      estimate = estimates,
      lower = lowers,
      upper = uppers
  )
  
  return(ci)
  
}