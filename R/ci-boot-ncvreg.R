#' Title
#'
#' @param boot 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
ci.boot.ncvreg <- function(boot, alpha = 0.2, quiet = FALSE) {
  
  all_draws <- boot[["draws"]]
    
  lowers <- apply(all_draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  uppers <- apply(all_draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
  center <- apply(all_draws, 2, function(x) quantile(x, 0.5, na.rm = TRUE))
  
  any_nas <- any(as.logical(apply(all_draws, 2, function(x) sum(is.na(x)) > 0)))
  if (any_nas & !quite) {warning("NAs in draws")}
  
  ci_info <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]), lower = lowers, upper = uppers, center = center, ci_method = ci_method)    
  
  return(ci_info)
 
}
