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
#' @param boot An object of type \code{boot_ncvreg}
#' @param alpha The desired significance level to obtain confidence intervals
#' with overall coverage greater than or equal to \code{1 - alpha}
#' @param quiet If \code{TRUE}, suppress warning that some bootstrap draws are
#' NA. This occurs when a given variable is \code{singlar} for a given bootstrap
#' sample (see \code{std}).
#'
#' @return An object with S3 class \code{data.frame} with the columns: \describe {
#' \item{variable}{If the original covariates \code{X} had column names,
#' these are included as the first column.}
#' \item{estimate}{A length \code{ncol(X)} vector of the estimates from the
#' lasso model fit on the original data corresponding to \code{lambda}.}
#' \item{lower}{A length \code{ncol(X)} vector of the lower bounds from the
#' Hybrid bootstrap draws with significance level \code{alpha}}
#' \item{upper}{A length \code{ncol(X)} vector of the upper bounds from the
#' Hybrid bootstrap draws with significance level \code{alpha}}}
#' 
#' @examples
#' 
#' data(Prostate)
#' boot <- boot_ncvreg(Prostate$X, Prostate$y)
#' confint(boot, alpha = 0.1)
#' 
#' @export
confint.boot_ncvreg <- function(boot, alpha = 0.2, quiet = FALSE) {
  
  if ((missing(boot) || class(boot) != "boot_ncvreg")) {
    stop("boot must be supplied and be of class boot_ncvreg.")
  }
  
  draws <- boot[["draws"]]
  
  any_nas <- any(as.logical(apply(draws, 2, function(x) sum(is.na(x)) > 0)))
  if (any_nas & !quiet) {
    warning("NAs in draws, this occurs when a varaibles is nonsingular (constant) for a given bootstrap draw.")
  }
  
  lowers <- apply(draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  uppers <- apply(draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
    
  ci <- data.frame(
      variable = names(boot[["estimates"]]),
      estimate = boot[["estimates"]],
      lower = lowers,
      upper = uppers
  )
  
  return(ci)
  
}