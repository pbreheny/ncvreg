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
#' @param boot An object of type \code{boot_ncvreg}
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
ci.boot_ncvreg <- function(boot, alpha = 0.2, quiet = FALSE, methods = "all") {
  
  method_list <- c("traditional", "posterior", "hybrid", "debiased", "debiased2")
  if (methods != "all") {
    method_list <- intersect(method_list, methods)
  }
  
  intervals_list <- list()
  
  if ("traditional" %in% method_list) {
    traditional_cis <- compute_intervals(boot[["point_estimates"]], alpha = alpha, quiet = quiet)
    intervals_list$traditional <- traditional_cis
  }
  
  if ("posterior" %in% method_list) {
    posterior_cis <- compute_intervals(boot[["fc_draws"]], alpha = alpha, quiet = quiet)
    intervals_list$posterior <- posterior_cis
  }
  
  if ("hybrid" %in% method_list) {
    point_estimates <- boot[["point_estimates"]]
    fc_draws <- boot[["fc_draws"]]
    
    result <- point_estimates
    result[is.na(point_estimates)] <- NA
    result[!is.na(point_estimates) & point_estimates == 0] <- fc_draws[!is.na(point_estimates) & point_estimates == 0]
    result[!is.na(point_estimates) & point_estimates != 0] <- point_estimates[!is.na(point_estimates) & point_estimates != 0]

    hybrid_cis <- compute_intervals(result, alpha = alpha, quiet = quiet)
    intervals_list$hybrid <- hybrid_cis
  }
  
  if ("debiased" %in% method_list) {
    debiased_cis <- compute_intervals(boot[["partial_correlations"]], alpha = alpha, quiet = quiet)
    intervals_list$debiased <- debiased_cis
  }
  
  if ("debiased2" %in% method_list) {
    point_estimates <- boot[["partial_correlations"]]
    debiased_draws <- boot[["debiased_draws"]]
    
    result <- point_estimates
    result[is.na(point_estimates)] <- NA
    result[!is.na(point_estimates) & point_estimates == 0] <- debiased_draws[!is.na(point_estimates) & point_estimates == 0]
    result[!is.na(point_estimates) & point_estimates != 0] <- point_estimates[!is.na(point_estimates) & point_estimates != 0]
    
    debiased2_cis <- compute_intervals(result, alpha = alpha, quiet = quiet)
    intervals_list$debiased2 <- debiased2_cis
  }
  
  ci_info_all <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]))
  ci_info <- list()
  
  for (method in names(intervals_list)) {
    
    ci_info[[method]] <- data.frame(
        lower = intervals_list[[method]][[1]],
        upper = intervals_list[[method]][[2]],
        variable = names(boot[["estimates"]])
      ) %>%
      mutate(method = method)
    
  }
  
  ci_info_all <- left_join(ci_info_all, do.call(rbind, ci_info), by = "variable")
  
  return(ci_info_all)
  
}
compute_intervals <- function(draws, alpha = 0.2, quiet = FALSE) {
  
  any_nas <- any(as.logical(apply(draws, 2, function(x) sum(is.na(x)) > 0)))
  if (any_nas & !quiet) {
    warning("NAs in draws")
  }
  
  lowers <- apply(draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  uppers <- apply(draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
  
  return(list(lowers, uppers))
  
}
