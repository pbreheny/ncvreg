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
  
  if (method == "quantile") {
    all_draws <- eb_boot[["draws"]]
    
    lowers <- apply(all_draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
    uppers <- apply(all_draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
    
    ci_info <- data.frame(estimate = eb_boot[["estimates"]], variable = names(eb_boot[["estimates"]]), lower = lowers, upper = uppers, method = "Lasso Boot")    
  } else if (method == "bca") {
    ci_info <- data.frame(
      "estimate" = eb_boot$estimates,
      "variable" = names(eb_boot$estimates),
      "ci" = t(sapply(1:ncol(eb_boot$draws), function(x) BCa_ci(eb_boot$draws[,x], eb_boot$estimates[x], alpha = alpha))),
      "method" = "Lasso Boot"
    )
    
    colnames(ci_info)[3:4] <- c("lower", "upper")
  }

  
  return(ci_info)
  
}

BCa_ci <- function(bootstrap_samples, original_estimate, alpha = 0.05) {
  # Calculate the number of bootstrap samples
  n <- length(bootstrap_samples)
  
  # Calculate bias correction (z0)
  z0 <- qnorm(mean(bootstrap_samples < original_estimate))
  
  # Compute acceleration (acc)
  j <- jackknife(bootstrap_samples, original_estimate)
  acc <- sum((mean(j) - j)^3) / (6 * sum((mean(j) - j)^2)^1.5)
  
  # Adjust alpha for two-tailed test
  alpha1 <- alpha / 2
  alpha2 <- 1 - alpha / 2
  
  # Compute adjusted percentiles
  p1 <- pnorm(z0 + (z0 + qnorm(alpha1)) / (1 - acc * (z0 + qnorm(alpha1))))
  p2 <- pnorm(z0 + (z0 + qnorm(alpha2)) / (1 - acc * (z0 + qnorm(alpha2))))
  
  # Calculate BCa confidence interval
  bca_ci <- quantile(bootstrap_samples, c(p1, p2))
  
  return(bca_ci)
}

jackknife <- function(bootstrap_samples, original_estimate) {
  # Compute jackknife estimates
  n <- length(bootstrap_samples)
  jk_estimates <- numeric(n)
  for (i in 1:n) {
    jk_sample <- bootstrap_samples[-i]
    jk_estimates[i] <- mean(jk_sample)
  }
  return(jk_estimates)
}