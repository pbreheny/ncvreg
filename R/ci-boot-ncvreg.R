#' Title
#'
#' @param eb_boot 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
ci.boot.ncvreg <- function(eb_boot, quiet = FALSE, method = "quantile", alpha = 0.2, original_data = NULL) {
  
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
  } else if (method == "bca") {
    ci_info <- data.frame(
      estimate = eb_boot[["estimates"]],
      variable = names(eb_boot[["estimates"]]),
      BCa_ci(eb_boot$draws, eb_boot$estimates, original_data = original_data, lambda = eb_boot$lambda, sigma2 = eb_boot$sigma2, alpha = alpha),
      method = "Lasso Boot")
    
    colnames(ci_info)[3:4] <- c("lower", "upper")
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
BCa_ci <- function(bootstrap_samples, original_estimate, original_data, lambda, sigma2, alpha = 0.05) {
  
  # Calculate the number of bootstrap samples
  n <- nrow(bootstrap_samples)
  
  # Calculate bias correction (z0)
  original_estimate <- bootf(original_data$X, original_data$y, lambda = lambda, sigma2 = sigma2, resample = FALSE, quantiles = "debiased")$draws
  quant <- sapply(1:ncol(bootstrap_samples), function(x) mean(bootstrap_samples[,x] < original_estimate[x]))
  quant <- ifelse(quant == 1, 1 - (1/n), quant)
  quant <- ifelse(quant == 0, 1/n, quant)
  z0 <- qnorm(quant)
  
  # acc
  acc <- jacknife_acc(original_data, lambda, sigma2)
  
  # Adjust alpha for two-tailed test
  alpha1 <- alpha / 2
  alpha2 <- 1 - alpha / 2
  
  # Compute adjusted percentiles
  p1 <- pnorm(z0 + (z0 + qnorm(alpha1)) / (1 - acc * (z0 + qnorm(alpha1))))
  p2 <- pnorm(z0 + (z0 + qnorm(alpha2)) / (1 - acc * (z0 + qnorm(alpha2))))
  
  # Calculate BCa confidence interval
  bca_ci_lower <- sapply(1:ncol(bootstrap_samples), function(x) quantile(bootstrap_samples[,x], p1[x]))
  bca_ci_upper <- sapply(1:ncol(bootstrap_samples), function(x) quantile(bootstrap_samples[,x], p2[x]))
  bca_ci <- cbind(bca_ci_lower, bca_ci_upper)
  
  return(bca_ci)
}
jacknife_acc <- function(original_data, lambda, sigma2) {
  
  # Compute jackknife estimates
  n <- nrow(original_data$X)
  jk_estimates <- numeric(n)
  
  # jk samples
  jk_samples <- matrix(nrow = n, ncol = ncol(original_data$X))
  for (i in 1:n) {
    jk_samples[i,] <- bootf_nosample(original_data$X[-i,], original_data$y[-i], lambda = lambda, sigma2 = sigma2, quantiles = "debiased")$draws
  }
  
  acc <- apply(jk_samples, 2, acceleration)
  
  return(acc)
}
acceleration <- function(j) {
  return(sum((mean(j) - j)^3) / (6 * sum((mean(j) - j)^2)^1.5))
}