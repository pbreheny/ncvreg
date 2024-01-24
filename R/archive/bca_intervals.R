bootf_nosample <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, quantiles = "disturbed") {

  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }

  p <- ncol(XX)
  n <- length(y)

  modes <- numeric(p)

  ynew <- y
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX)
  nonsingular <- attr(xnew, "nonsingular")

  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  full_rescale_factor <- rescale * rescaleX

  lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
  lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
  if (lambda_min > lambda_max | lambda > lambda_max) {
    lambda_max <- lambda + lambda / 100
    nlambda <- 2
  }
  ## Could use better logic to speed up
  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
  ncvreg.args$X <- xnew
  ncvreg.args$y <- ynew
  ncvreg.args$penalty <- "lasso"
  ncvreg.args$lambda <- lambda_seq

  ## Ignores user specified lambda.min and nlambda
  # fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])
  fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)

  coefs <- coef(fit, lambda = lambda)
  modes <- coefs[-1] ## Coefs only returned for nonsingular columns of X
  partial_residuals <-  ynew - (coefs[1] + as.numeric(xnew %*% modes) - (xnew * matrix(modes, nrow=nrow(xnew), ncol=ncol(xnew), byrow=TRUE)))

  z <- (1/n)*colSums(xnew * partial_residuals)
  se <- sqrt(sigma2 / n)

  ## Tails I am transferring on to (log probability in each tail)
  obs_lw <- pnorm(0, z + lambda, se, log.p = TRUE)
  obs_up <- pnorm(0, z - lambda, se, lower.tail = FALSE, log.p = TRUE)

  ## alt
  obs_p_lw <- obs_lw + ((z*lambda*n) / sigma2)
  obs_p_up <- obs_up - ((z*lambda*n) / sigma2)

  ## But I need to use this to find the proportion of each to the overall probability
  frac_lw_log <- ifelse(is.infinite(exp(obs_p_lw - obs_p_up)), 0, obs_p_lw - obs_p_up - log(1 + exp(obs_p_lw - obs_p_up)))
  frac_up_log <- ifelse(is.infinite(exp(obs_p_up - obs_p_lw)), 0, obs_p_up - obs_p_lw - log(1 + exp(obs_p_up - obs_p_lw)))

  dmodes <- ifelse(
    modes <= 0,
    dnorm(modes, z + lambda, se, log = TRUE) + obs_lw - frac_lw_log,
    dnorm(modes, z - lambda, se, log = TRUE) + obs_up - frac_up_log
  )


  spans <- runif(length(modes), 0, 3)

  accepted <- logical(length(modes))
  draws <- matrix(ncol = p, nrow = 1)
  iters <- 0
  while (any(!accepted) & iters < 1000) {
    for (i in 1:length(dmodes)) {
      if (!accepted[i]) {
        curr_sign <- sample(c(-1, 1), 1)
        curr_x <- modes[i] + curr_sign * spans[i] * se
        curr_thresh <- log(runif(1))


        curr_dens <- ifelse(
          curr_x <= 0,
          dnorm(curr_x, z[i] + lambda, se, log = TRUE) + obs_lw[i] - frac_lw_log[i] - dmodes[i],
          dnorm(curr_x, z[i] - lambda, se, log = TRUE) + obs_up[i] - frac_up_log[i] - dmodes[i]
        )

        if (curr_dens >= curr_thresh) {
          draws[1,i] <- curr_x
          accepted[i] <- TRUE
        } else if (sign(curr_x) != sign(modes[i])) {
          curr_x <- modes[i] + curr_sign * -1 * spans[i] * se
          curr_dens <- ifelse(
            curr_x <= 0,
            dnorm(curr_x, z[i] + lambda, se, log = TRUE) + obs_lw[i] - frac_lw_log[i] - dmodes[i],
            dnorm(curr_x, z[i] - lambda, se, log = TRUE) + obs_up[i] - frac_up_log[i] - dmodes[i]
          )
          if (curr_dens >= curr_thresh) {
            draws[1,i] <- curr_x
            accepted[i] <- TRUE
          }
        }

        spans[i] <- runif(1, 0, spans[i])
      }
    }
    iters <- iters + 1
  }

  draws <- draws[,nonsingular,drop=FALSE] * full_rescale_factor

  modes[nonsingular] <- (modes * rescale) * rescaleX
  if (length(nonsingular) < ncol(draws)) draws[,!(1:ncol(draws) %in% nonsingular)] <- NA
  modes[!(1:length(modes) %in% nonsingular)] <- NA

  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  return(ret)

}
bootf_nosample2 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, quantiles = "zs") {

  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }

  p <- ncol(XX)
  n <- length(y)

  modes <- numeric(p)

  ynew <- y
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX)
  nonsingular <- attr(xnew, "nonsingular")

  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  full_rescale_factor <- rescale * rescaleX

  lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
  lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
  if (lambda_min > lambda_max | lambda > lambda_max) {
    lambda_max <- lambda + lambda / 100
    nlambda <- 2
  }
  ## Could use better logic to speed up
  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
  ncvreg.args$X <- xnew
  ncvreg.args$y <- ynew
  ncvreg.args$penalty <- "lasso"
  ncvreg.args$lambda <- lambda_seq

  ## Ignores user specified lambda.min and nlambda
  # fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])
  fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)

  coefs <- coef(fit, lambda = lambda)

  modes <- coefs[-1] ## Coefs only returned for nonsingular columns of X
  partial_residuals <-  ynew - (coefs[1] + as.numeric(xnew %*% modes) - (xnew * matrix(modes, nrow=nrow(xnew), ncol=ncol(xnew), byrow=TRUE)))

  z <- (1/n)*colSums(xnew * partial_residuals)
  draws <- matrix(ncol = p, nrow = 1)
  draws[1,] <- z[nonsingular] * full_rescale_factor

  modes[nonsingular] <- (modes * rescale) * rescaleX

  if (length(nonsingular) < ncol(draws)) draws[,!(1:ncol(draws) %in% nonsingular)] <- NA
  modes[!(1:length(modes) %in% nonsingular)] <- NA

  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")

  return(ret)

}
bootf_nosample3 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, quantiles = "mode") {

  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }

  p <- ncol(XX)
  n <- length(y)

  modes <- numeric(p)
  ynew <- y
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX)
  nonsingular <- attr(xnew, "nonsingular")

  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  full_rescale_factor <- rescale * rescaleX

  lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
  lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
  if (lambda_min > lambda_max | lambda > lambda_max) {
    lambda_max <- lambda + lambda / 100
    nlambda <- 2
  }
  ## Could use better logic to speed up
  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
  ncvreg.args$X <- xnew
  ncvreg.args$y <- ynew
  ncvreg.args$penalty <- "lasso"
  ncvreg.args$lambda <- lambda_seq

  ## Ignores user specified lambda.min and nlambda
  fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])
  # fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)

  coefs <- coef(fit, lambda = lambda)

  modes <- coefs[-1] ## Coefs only returned for nonsingular columns of X
  modes[nonsingular] <- (modes * rescale) * rescaleX
  modes[!(1:length(modes) %in% nonsingular)] <- NA

  ret <- list(modes, modes)
  names(ret) <- c("draws", "modes")

  return(ret)

}
bootf_nosample5 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, quantiles = "zerosample") {

  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }

  p <- ncol(XX)
  n <- length(y)

  ynew <- y
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX)
  nonsingular <- attr(xnew, "nonsingular")

  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1) ## may need to take a closer look at this
  } else {
    rescaleX <- 1
  }
  full_rescale_factor <- rescale * rescaleX


  lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
  lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
  if (lambda_min > lambda_max | lambda > lambda_max) {
    lambda_max <- lambda + lambda / 100
    nlambda <- 2
  }

  ## Could use better logic to speed up
  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
  ncvreg.args$X <- xnew
  ncvreg.args$y <- ynew
  ncvreg.args$penalty <- "lasso"
  ncvreg.args$lambda <- lambda_seq

  ## Ignores user specified lambda.min and nlambda
  fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])
  # fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)

  coefs <- coef(fit, lambda = lambda)
  modes <- coefs[-1] ## Coefs only returned for nonsingular columns of X
  partial_residuals <-  ynew - (coefs[1] + as.numeric(xnew %*% modes) - (xnew * matrix(modes, nrow=nrow(xnew), ncol=ncol(xnew), byrow=TRUE)))

  z <- (1/n)*colSums(xnew * partial_residuals)
  se <- sqrt(sigma2 / n)

  ## Tails I am transferring on to (log probability in each tail)
  obs_lw <- pnorm(0, z + lambda, se, log.p = TRUE)
  obs_up <- pnorm(0, z - lambda, se, lower.tail = FALSE, log.p = TRUE)

  ## alt
  obs_p_lw <- obs_lw + ((z*lambda*n) / sigma2)
  obs_p_up <- obs_up - ((z*lambda*n) / sigma2)

  ## But I need to use this to find the proportion of each to the overall probability
  frac_lw_log <- ifelse(is.infinite(exp(obs_p_lw - obs_p_up)), 0, obs_p_lw - obs_p_up - log(1 + exp(obs_p_lw - obs_p_up)))
  frac_up_log <- ifelse(is.infinite(exp(obs_p_up - obs_p_lw)), 0, obs_p_up - obs_p_lw - log(1 + exp(obs_p_up - obs_p_lw)))

  ## Can make much more efficient
  ps <- runif(length(frac_lw_log))
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  draws <- matrix(ncol = p, nrow = 1)
  draws[1,nonsingular] <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(log_one_minus_ps + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  ) * full_rescale_factor
  draws[1, nonsingular[modes != 0]] <- modes[modes != 0] * full_rescale_factor[modes != 0]


  tmp <- modes
  modes <- numeric(p)
  modes[nonsingular] <- tmp * full_rescale_factor
  if (length(nonsingular) < ncol(draws)) draws[1,!(1:ncol(draws) %in% nonsingular)] <- NA
  modes[!(1:length(modes) %in% nonsingular)] <- NA

  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  return(ret)

}
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
  } else if (method == "bca2") {
    ci_info <- data.frame(
      "estimate" = eb_boot$estimates,
      "variable" = names(eb_boot$estimates),
      "ci" = BCa_ci2(eb_boot$draws, eb_boot$estimates, original_data = original_data, lambda = eb_boot$lambda, sigma2 = eb_boot$sigma2, alpha = alpha),
      "method" = "Lasso Boot"
    )
    
    colnames(ci_info)[3:4] <- c("lower", "upper")
  } else if (method == "bca5") {
    ci_info <- data.frame(
      "estimate" = eb_boot$estimates,
      "variable" = names(eb_boot$estimates),
      "ci" = BCa_ci5(eb_boot$draws, eb_boot$estimates, original_data = original_data, lambda = eb_boot$lambda, sigma2 = eb_boot$sigma2, alpha = alpha),
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
  
  # Calculate BCa confidence interval
  bca_ci_lower <- sapply(1:ncol(bootstrap_samples), function(x) quantile(bootstrap_samples[,x], p1[x]))
  bca_ci_upper <- sapply(1:ncol(bootstrap_samples), function(x) quantile(bootstrap_samples[,x], p2[x]))
  bca_ci <- cbind(bca_ci_lower, bca_ci_upper)
  
  return(bca_ci)
}

BCa_ci2 <- function(bootstrap_samples, original_estimate, original_data, lambda, sigma2, alpha = 0.05) {
  
  # Calculate the number of bootstrap samples
  n <- nrow(bootstrap_samples)
  
  # Calculate bias correction (z0)
  original_estimate <- bootf_nosample2(original_data$X, original_data$y, lambda = lambda, sigma2 = sigma2)$draws
  quant <- sapply(1:ncol(bootstrap_samples), function(x) mean(bootstrap_samples[,x] < original_estimate[x]))
  quant <- ifelse(quant == 1, 1 - (1/n), quant)
  quant <- ifelse(quant == 0, 1/n, quant)
  z0 <- qnorm(quant)
  
  # acc
  acc <- jacknife_acc2(original_data, lambda, sigma2)
  
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

BCa_ci5 <- function(bootstrap_samples, original_estimate, original_data, lambda, sigma2, alpha = 0.05) {
  
  # Calculate the number of bootstrap samples
  n <- nrow(bootstrap_samples)
  
  # Calculate bias correction (z0)
  quant <- sapply(1:ncol(bootstrap_samples), function(x) mean(bootstrap_samples[,x] < original_estimate[x]))
  quant <- ifelse(quant == 1, 1 - (1/n), quant)
  quant <- ifelse(quant == 0, 1/n, quant)
  z0 <- qnorm(quant)
  
  # acc
  acc <- jacknife_acc5(original_data, lambda, sigma2)
  
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
    jk_samples[i,] <- bootf_nosample(original_data$X[-i,], original_data$y[-i], lambda = lambda, sigma2 = sigma2)$draws
  }
  
  acc <- apply(jk_samples, 2, acceleration)
  
  return(acc)
}

jacknife_acc2 <- function(original_data, lambda, sigma2) {
  
  # Compute jackknife estimates
  n <- nrow(original_data$X)
  jk_estimates <- numeric(n)
  
  # jk samples
  jk_samples <- matrix(nrow = n, ncol = ncol(original_data$X))
  for (i in 1:n) {
    jk_samples[i,] <- bootf_nosample2(original_data$X[-i,], original_data$y[-i], lambda = lambda, sigma2 = sigma2)$draws
  }
  
  acc <- apply(jk_samples, 2, acceleration)
  
  return(acc)
}

jacknife_acc5 <- function(original_data, lambda, sigma2) {
  
  # Compute jackknife estimates
  n <- nrow(original_data$X)
  jk_estimates <- numeric(n)
  
  # jk samples
  jk_samples <- matrix(nrow = n, ncol = ncol(original_data$X))
  for (i in 1:n) {
    jk_samples[i,] <- bootf_nosample5(original_data$X[-i,], original_data$y[-i], lambda = lambda, sigma2 = sigma2)$draws
  }
  
  acc <- apply(jk_samples, 2, acceleration)
  
  return(acc)
}
