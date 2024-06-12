#' Calculating percentile confidence intervals from bootstrap draws
#' 
#' Calculates percentile based confidence intervals from hybrid bootstrap draws
#' obtained from \code{boot_ncvreg}. This produces intervals that are generally
#' conservative on average, especially when n is small or the underlying covarites are
#' especially sparse. However, intervals may undercover when correlation amoung
#' the predictors is especially high. Although overall coverage is generally
#' the intervals produced do no attempt to debias and as a result generally
#' undercover large covariates and overcover covariates near zero.
#'
#' @param bootfit An object of type \code{boot_ncvreg}
#' @param alpha The desired significance level to obtain confidence intervals
#' with overall coverage greater than or equal \code{1 - alpha}
#' @param quiet If \code{TRUE}, suppress warning that some bootstrap draws are
#' NA. This occurs when a given variable is \code{singlar} for a given bootstrap
#' sample (see \code{std}).
#'
#' @return An object with S3 class \code{data.frame} with the columns:
#' \describe { \item{estimate}{A length \code{ncol(X)} vector of the estimates from the
#' lasso model fit on the original data corresponding to \code{lambda}.}}
#' @export
#'
#' @examples
ci.boot_ncvreg <- function(boot, alpha = 0.2, quiet = FALSE) {
  
  all_draws <- boot[["draws"]]
    
  lowers <- apply(all_draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  uppers <- apply(all_draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
  center <- apply(all_draws, 2, function(x) quantile(x, 0.5, na.rm = TRUE))
  
  any_nas <- any(as.logical(apply(all_draws, 2, function(x) sum(is.na(x)) > 0)))
  if (any_nas & !quiet) {warning("NAs in draws")}
  
  ci_info <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]), lower = lowers, upper = uppers, center = center)    
  
  return(ci_info)
 
}
