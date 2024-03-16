#' Title
#'
#' @param boot 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
ci.boot.ncvreg <- function(boot, quiet = FALSE, ci_method = "quantile", alpha = 0.2, original_data = NULL) {
  
  all_draws <- boot[["draws"]]
  if (ci_method == "quantile") {
    
    lowers <- apply(all_draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
    uppers <- apply(all_draws, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
    
    ci_info <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]), lower = lowers, upper = uppers, ci_method = ci_method)    
  
  } else if (ci_method == "bucketfill") {
  
    estimates <- boot[["estimates"]]
    bounds <- do.call(rbind, lapply(1:ncol(all_draws), function(x) fill_bucket(all_draws[,x], estimates[x], alpha)))
    ci_info <- data.frame(estimate = estimates, variable = names(boot[["estimates"]]), lower = bounds[,1], upper = bounds[,2], ci_method = ci_method)    
  
  } else if (ci_method == "identity") {
    
    lowers <- all_draws[1,]
    uppers <- all_draws[2,]
    
    ci_info <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]), lower = lowers, upper = uppers, ci_method = ci_method)    
  } else if (ci_method == "bca") {
    ci_info <- data.frame(
      estimate = boot[["estimates"]],
      variable = names(boot[["estimates"]]),
      BCa_ci(all_draws, original_data = original_data, lambda = boot$lambda, sigma2 = boot$sigma2, alpha = alpha, method = boot$method),
      ci_method = ci_method)
    
    colnames(ci_info)[3:4] <- c("lower", "upper")
  } else if (ci_method == "mvn_uni") {
    
    means <- apply(all_draws, 2, mean)
    vars <- (colSums(scale(all_draws, scale = FALSE)^2) / (nrow(all_draws) - 1)) 
    
    # tmp <- mapply(draw_samples, mean=means, variance=vars, MoreArgs=list(n=10000))
    # cis <- apply(tmp, 2, function(x) quantile(x, c(alpha / 2, 1 - (alpha/2))))
    df <- max(mean(apply(boot$modes, 1, function(x) (sum(x != 0)))), 1)
    cis <- mapply(produce_t_cis, mean=means, variance=vars, MoreArgs=list(alpha=alpha, df = df))
    
    ci_info <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]), lower = cis[1,], upper = cis[2,], ci_method = ci_method)  
    
  } else if (ci_method == "mvn_corrected") {
    
    means <- apply(all_draws, 2, mean)
    vars <- colSums(scale(all_draws, scale = FALSE)^2) / (nrow(all_draws) - 1)
    
    rate <- (nrow(original_data$X) * boot$lambda) / boot$sigma2
    # rate <- boot$lambda
    rescale <- attr(ncvreg::std(original_data$X), "scale") ## need to think about this more
    tmp <- mapply(draw_samples_corrected, mean=means, variance=vars, rescale=rescale, MoreArgs=list(n=10000, rate = rate))
    # tmp <- mapply(draw_samples_corrected, idx=1:ncol(all_draws), rescale=rescale, MoreArgs=list(rate = rate, all_draws = all_draws))
   
    cis <- apply(tmp, 2, function(x) quantile(x, c(alpha / 2, 1 - (alpha/2))))
    ci_info <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]), lower = cis[1,], upper = cis[2,], ci_method = ci_method)    
  } else if (ci_method == "full_debias") {
    
    rescale <- attr(ncvreg::std(original_data$X), "scale")^(-1) ## need to think about this more
    cis <- sapply(1:ncol(original_data$X), full_debias, lassoboot = boot, original_data = original_data, alpha = alpha, rescale = rescale)
    ci_info <- data.frame(estimate = boot[["estimates"]], variable = names(boot[["estimates"]]), lower = cis[1,], upper = cis[2,], ci_method = ci_method)    
    
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
BCa_ci <- function(bootstrap_samples, original_data, lambda, sigma2, alpha = 0.05, method) {
  
  # Calculate the number of bootstrap samples
  n <- nrow(bootstrap_samples)
  
  # Calculate bias correction (z0)
  original_estimate <- bootf(original_data$X, original_data$y, lambda = lambda, sigma2 = sigma2, resample = FALSE, method = method)$draws
  quant <- sapply(1:ncol(bootstrap_samples), function(x) mean(bootstrap_samples[,x] < original_estimate[x]))
  quant <- ifelse(quant == 1, 1 - (1/n), quant)
  quant <- ifelse(quant == 0, 1/n, quant)
  z0 <- qnorm(quant)
  
  # acc
  acc <- jacknife_acc(original_data, lambda, sigma2, method = method)
  
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
jacknife_acc <- function(original_data, lambda, sigma2, method) {
  
  # Compute jackknife estimates
  n <- nrow(original_data$X)
  jk_estimates <- numeric(n)
  
  # jk samples
  jk_samples <- matrix(nrow = n, ncol = ncol(original_data$X))
  for (i in 1:n) {
    jk_samples[i,] <- bootf(original_data$X[-i,], original_data$y[-i], lambda = lambda, sigma2 = sigma2, method = method, resample = FALSE)$draws
  }
  
  acc <- apply(jk_samples, 2, acceleration)
  
  return(acc)
}
acceleration <- function(j) {
  return(sum((mean(j) - j)^3) / (6 * sum((mean(j) - j)^2)^1.5))
}
draw_samples <- function(mean, variance, n) {
  rnorm(n, mean, sqrt(variance))
}
produce_normal_cis <- function(mean, variance, alpha) {
  mean + c(-1, 1)*sqrt(variance)*qnorm(1-alpha/2)
}
produce_t_cis <- function(mean, variance, alpha, df) {
  mean + c(-1, 1)*sqrt(variance)*qt(1-alpha/2, df = df)
}
## Make sure everything is on agreeable scale
draw_samples_corrected <- function(mean, variance, rescale, n, rate) {
  tmp <- rnorm(n, mean, sqrt(variance))
  probs <- (rate / 2) * exp(-rate*abs(tmp*rescale))
  probs <- probs / sum(probs)
  sample(tmp, replace = TRUE, prob = probs)
}
# draw_samples_corrected <- function(idx, rescale, rate, all_draws) {
#   probs <- (rate / 2) * exp(-rate*abs(all_draws[,idx]*rescale))
#   probs <- probs / sum(probs)
#   sample(all_draws[,idx], replace = TRUE, prob = probs)
# }
full_debias <- function(which_var, lassoboot, original_data, alpha, rescale) {
  ds <- lassoboot$draws[,which_var]
  ms <- lassoboot$modes[,which_var]
  
  
  
  lam <- lassoboot$lambda * rescale
  sdv <- lassoboot$sigma2
  
  # mns <- ifelse(ms == 0, ds, lam) * sign(ds)
  # corrections <- rnorm(n = length(mns), mns, sd = sqrt(sdv / nrow(original_data$X)))
  corrections <- mns
  return(quantile(ds + corrections, c(alpha / 2, 1 - (alpha / 2))))
}
