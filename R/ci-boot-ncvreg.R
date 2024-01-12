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
  
  if (method == "quantile") {
    all_draws <- eb_boot[["draws"]]
    
    lowers <- apply(all_draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
    uppers <- apply(all_draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
    
    ci_info <- data.frame(estimate = eb_boot[["estimates"]], variable = names(eb_boot[["estimates"]]), lower = lowers, upper = uppers, method = "Lasso Boot")    
  } else if (method == "bca") {
    ci_info <- data.frame(
      "estimate" = eb_boot$estimates,
      "variable" = names(eb_boot$estimates),
      "ci" = BCa_ci(eb_boot$draws, eb_boot$estimates, original_data = original_data, lambda = eb_boot$lambda, sigma2 = eb_boot$sigma2, alpha = alpha),
      "method" = "Lasso Boot"
    )
    
    colnames(ci_info)[3:4] <- c("lower", "upper")
  }

  
  return(ci_info)
  
}

BCa_ci <- function(bootstrap_samples, original_estimate, original_data, lambda, sigma2, alpha = 0.05) {
  
  # Calculate the number of bootstrap samples
  n <- nrow(bootstrap_samples)
  
  # Calculate bias correction (z0)
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
  print(p1)
  print(p2)
  
  # Calculate BCa confidence interval
  bca_ci_lower <- sapply(1:ncol(bootstrap_samples), function(x) quantile(bootstrap_samples[,x], p1[x]))
  bca_ci_upper <- sapply(1:ncol(bootstrap_samples), function(x) quantile(bootstrap_samples[,x], p2[x]))
  bca_ci <- cbind(bca_ci_lower, bca_ci_upper)
  print(bca_ci)
  
  return(bca_ci)
}

acceleration <- function(j) {
  return(sum((mean(j) - j)^3) / (6 * sum((mean(j) - j)^2)^1.5))
}

jacknife_acc <- function(original_data, lambda, sigma2) {
  
  # Compute jackknife estimates
  n <- nrow(original_data$X)
  jk_estimates <- numeric(n)
  
  # jk samples
  jk_samples <- matrix(nrow = n, ncol = ncol(original_data$X))
  for (i in 1:n) {
    jk_samples[i,] <- bootf_nosample(original_data$X, original_data$y, lambda = lambda, sigma2 = sigma2)$draws
  }
  
  acc <- apply(jk_samples, 2, acceleration)
  
  return(acc)
}
