#' Title
#'
#' @param eb_boot 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
ci.boot.ncvreg.r <- function(eb_boot, quiet = FALSE) {
  
  rm_draw <- apply(eb_boot[["draws"]], 2, function(x) sum(is.na(x))); names(rm_draw) <- names(eb_boot$estimates)
  rm_draw <- rm_draw[rm_draw != 0]
  
  # all_draws <- rbind(eb_boot[["lowers"]], eb_boot[["uppers"]])
  all_draws <- eb_boot[["draws"]]
  # lowers <- apply(eb_boot[["lowers"]], 2, mean, na.rm = TRUE)
  # uppers <- apply(eb_boot[["uppers"]], 2, mean, na.rm = TRUE)
  lowers <- apply(all_draws, 2, function(x) quantile(x, 0.1, na.rm = TRUE))
  uppers <- apply(all_draws, 2, function(x) quantile(x, 0.9, na.rm = TRUE))
  
  ci_info <- data.frame(estimate = eb_boot[["estimates"]], variable = names(eb_boot[["estimates"]]), lower = lowers, upper = uppers, method = "Lasso Boot")
  
  return(ci_info)
  
}