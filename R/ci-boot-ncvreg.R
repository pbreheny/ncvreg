#' Title
#'
#' @param eb_boot 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
ci.boot.ncvreg <- function(eb_boot, quiet = FALSE, method = "quantile", alpha = 0.2) {
  
  all_draws <- eb_boot[["draws"]]
  if (method == "quantile") {
    
    lowers <- apply(all_draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
    uppers <- apply(all_draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
    
    ci_info <- data.frame(estimate = eb_boot[["estimates"]], variable = names(eb_boot[["estimates"]]), lower = lowers, upper = uppers, method = "Lasso Boot")    
  
  } else if (method == "bucketfill") {
  
    estimates <- eb_boot[["estimates"]]
    bounds <- do.call(rbind, lapply(1:ncol(all_draws), function(x) fill_bucket(all_draws[,x], estimates[x], alpha)))
    ci_info <- data.frame(estimate = estimates, variable = names(eb_boot[["estimates"]]), lower = bounds[,1], upper = bounds[,2], method = "Lasso Boot")    
  
  } else if (method == "identity") {
    
    lowers <- all_draws[1,]
    uppers <- all_draws[2,]
    
    ci_info <- data.frame(estimate = eb_boot[["estimates"]], variable = names(eb_boot[["estimates"]]), lower = lowers, upper = uppers, method = "Lasso Boot")    
  }
  
  return(ci_info)
  
}
fill_bucket <- function(draws, estimate, alpha) {
  
  draws <- na.omit(draws)
  nd <- length(draws)
  ee <- sum(draws == estimate)
  to_allocate <- ceiling(nd*(1-alpha)) - ee
  # if (to_allocate <= 0) stop("No draws to allocate, too many draws equal estimate")
  
  to_allocate_below <- floor(to_allocate / 2) + (estimate < 0)*(to_allocate %% 2)
  to_allocate_above <- floor(to_allocate / 2) + (estimate > 0)*(to_allocate %% 2)
  ## losing 1 if equal zero... need to think about this more
  
  add_lower <- max(c(to_allocate_above - sum(draws > estimate), 0))
  add_upper <- max(c(to_allocate_below - sum(draws < estimate), 0))
  
  to_allocate_below <- to_allocate_below + add_lower
  to_allocate_above <- to_allocate_above + add_upper
  
  lower <- ifelse(sum(draws < estimate) >= to_allocate_below, 1 + sum(draws < estimate) - to_allocate_below, 1)
  upper <- ifelse(sum(draws > estimate) >= to_allocate_above, nd - sum(draws > estimate) + to_allocate_above, nd)
  
  bounds <- sort(draws)[c(lower, upper)]
  return(data.frame(lower = bounds[1], upper = bounds[2]))
  
} 
