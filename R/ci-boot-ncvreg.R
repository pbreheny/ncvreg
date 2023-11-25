#' Title
#'
#' @param eb_boot 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
ci.boot.ncvreg <- function(eb_boot, quiet = FALSE) {
  
  rm_lower <- apply(eb_boot[["lowers"]], 2, function(x) sum(is.na(x))); names(rm_lower) <- names(eb_boot$estimates)
  rm_upper <- apply(eb_boot[["uppers"]], 2, function(x) sum(is.na(x))); names(rm_upper) <- names(eb_boot$estimates)
  
  rm_lower <- rm_lower[rm_lower != 0]
  rm_upper <- rm_upper[rm_upper != 0]
  
  if (length(rm_lower) > 0 ) {
    if (!quiet) {
      print(paste0(length(rm_lower), " total variables with NA entries, summary: "))
      print(summary(rm_lower))
    } 
  }
  
  # all_draws <- rbind(eb_boot[["lowers"]], eb_boot[["uppers"]])
  lowers <- apply(eb_boot[["lowers"]], 2, mean, na.rm = TRUE)
  uppers <- apply(eb_boot[["uppers"]], 2, mean, na.rm = TRUE)
  # lowers <- apply(all_draws, 2, function(x) quantile(x, 0.1, na.rm = TRUE))
  # uppers <- apply(all_draws, 2, function(x) quantile(x, 0.9, na.rm = TRUE))
  
  ci_info <- data.frame(estimate = eb_boot[["estimates"]], variable = names(eb_boot[["estimates"]]), lower = lowers, upper = uppers, method = "Lasso Boot")
  
  return(ci_info)
  
}