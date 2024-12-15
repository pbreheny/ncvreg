#' Title
#'
#' @param X
#' @param y
#' @param fit
#' @param lambda
#' @param sigma
#' @param family
#' @param penalty
#' @param gamma
#' @param alpha
#' @param level
#' @param posterior
#'
#' @return
#' @export
#'
#' @examples
pipe <- function(X, y, fit, lambda, sigma,
                 family = c("gaussian", "binomial", "poisson"),
                 penalty = c("MCP", "SCAD", "lasso"),
                 gamma = switch(penalty, SCAD = 3.7, 3),
                 alpha = 1, level = 0.95,
                 posterior = FALSE
) {
  
  ## If user wants more control, call ncvreg or cv.ncvreg directly
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  
  if (missing(fit) & (missing(X) | missing(y))) {
    stop("Must supply a fit of class ncvreg or cv.ncvreg or both X and y.")
  }
  
  if ((!missing(fit) && is.null(fit$X)) & missing(X)) {
    stop("This procedure requires X. Either supply X, or fit the model using the option 'returnX = TRUE'")
  }
  
  if ((!missing(fit) && !is.null(fit$X)) && (!missing(X) && !identical(fit$X, X))) {
    stop(glue("X supplied both in {class(fit)} and as an argument to pipe and they are not the same. It is unclear which should be used."))
  }
  
  if (!missing(fit) && (!missing(y) && !identical(fit$y, y))) {
    stop(glue("y supplied along with object of class {class(fit)} which also contains y and they are not the same. It is unclear which should be used."))
  }
  
  if (alpha < 1) {
    stop("alpha less than 1 not yet implimented.")
  }
  
  ## Get a fit in some way, after this, every version of supplied X and y is the same
  if (missing(fit)) {
    if (missing(lambda)) {
      cv_fit <- cv.ncvreg(X, y, family = family, penalty = penalty, gamma = gamma, alpha = alpha, returnX = TRUE)
      fit <- cv_fit$fit
      lambda <- cv_fit$lambda.min
    } else {
      fit <- ncvreg(X, y, family = family, penalty = penalty, gamma = gamma, alpha = alpha, returnX = TRUE)
    }
  }
  
  if (missing(lambda)) {
    message("No lambda provided, using cv.ncvreg to select the value of lambda that minimizes CVE.")
    cv_fit <- cv.ncvreg(X, y, family = family, penalty = penalty, gamma = gamma, alpha = alpha, returnX = TRUE)
    lambda <- cv_fit$lambda.min
  }
  
  ## Get standardized X
  if (is.null(fit$X)) {
    if (is.null(attr(X, "scale"))) {
      XX <- ncvreg::std(X)
    } else {
      XX <- X
    }
  } else {
    XX <- fit$X
  }
  rescale_factorX <- attr(XX, "scale")
  yy <- fit$y ## not centered??
  if (fit$family == "gaussian") yy <- yy - mean(yy) ## center y
  
  ## Check for nonsingular?? Can't do unless have original X
  ## likely not a huge deal for non-bootstrap scenario
  # nonsingular <- attr(XX, "nonsingular")
  # if (length(nonsingular) != ncol(XX)) {
  #   warning("Columns in X are singular, inference will not be provided for these features")
  # }
  
  p <- ncol(XX)
  n <- nrow(XX)
  
  # initialize
  beta_PIPE <- numeric(p)
  sigma_PIPE <- numeric(p)
  
  if (fit$family == "gaussian"){
    
    # Load model
    beta <- coef(fit, lambda = lambda)[-1]
    beta <- beta * rescale_factorX
    S <- beta != 0
    S_hat <- which(S)
    N_hat <- which(!S)
    
    # Compute sigma if not provided
    if (missing(sigma)) {
      sigma <- sqrt(crossprod(yy - XX %*% beta) / (n - sum(beta != 0)))
    }
    
    yhat <- as.numeric(XX %*% beta)
    partial_residuals <- (yy - yhat) + (XX * matrix(beta, nrow = nrow(XX), ncol = ncol(XX), byrow = TRUE))
    beta_PIPE <- (1/n) * colSums(XX * partial_residuals)
    
    ## Compute PIPE Variance
    ## For features in the estimated support
    for (i in S_hat){
      
      S_hat_i <- S
      S_hat_i[i] <- FALSE
      
      ## Compute variance
      Xsi <- XX[,S_hat_i,drop = FALSE]
      Qsi <- diag(nrow(XX)) - tcrossprod(qr.Q(qr(Xsi)))
      adjusted_n <- t(XX[,i,drop = FALSE]) %*% Qsi %*% XX[,i,drop = FALSE]
      sigma_PIPE[i] <- sqrt(sigma^2 / adjusted_n)
      
    }
    
    # For null features
    if (length(S_hat) > 0 & length(N_hat) > 0) {
      Xs <- XX[,S_hat]
      Qs <- diag(nrow(XX)) - tcrossprod(qr.Q(qr(Xs)))
    } else {
      Qs <- diag(nrow(XX))
    }
    
    for (i in N_hat) {
      
      adjusted_n <- t(XX[,i,drop = FALSE]) %*% Qs %*% XX[,i,drop = FALSE]
      sigma_PIPE[i] <- sqrt(sigma^2 / adjusted_n)
      
    }
    
  } else {
    
    ## Refit the model with standardized X -> I don't love this... but not sure of a better solution
    if (!all(abs(attr(XX, "scale") - 1) < 1e-8)) {
      fit <- ncvreg::ncvreg(
        XX, yy, family = fit$family, penalty = fit$penalty, gamma = fit$gamma, alpha = fit$alpha
      )
    }
    
    pii <- predict(fit, XX, type = "response", lambda = lambda)
    beta <- coef(fit, lambda = lambda)
    S <- beta[-1] != 0
    S_hat <- which(S)
    N_hat <- which(!S)
    
    if (fit$family == "binomial") {
      A <- pii * (1 - pii)
    } else if (fit$family == "poisson") {
      A <- pii
    }
    
    v <- yy - pii
    W <- diag(A)
    Y_pse <- diag(1/A) %*% v + predict(fit, XX, type = "link", lambda = lambda)
    
    # compute weight matrix
    sqrtW <- sqrt(W)
    yy_pse <- sqrtW %*% Y_pse
    X_int <- cbind(rep(1, nrow(XX)), XX)
    weights <- colSums((W %*% XX) * XX) / n
    
    # Compute pipe for features in the estimated support
    for (i in S_hat) {
      
      S_hat_i <- S
      S_hat_i[i] <- FALSE
      
      Xs <- X_int[, c(TRUE, S_hat_i), drop = FALSE]
      XXs <- sqrtW %*% Xs
      xx <- crossprod(sqrtW, XX[,i])
      
      Qs <- diag(nrow(XX)) - tcrossprod(qr.Q(qr(XXs)))
      beta_PIPE[i] <- t(xx) %*% (yy_pse - XXs %*% beta[c(TRUE, S_hat_i)]) / crossprod(xx)
      sigma_PIPE[i] <- sqrt(1 / (t(xx) %*% Qs %*% xx))
      
    }
    
    # Compute pipe for null features
    Xs <- X_int[,c(TRUE, S), drop = FALSE]
    XXs <- sqrtW %*% Xs
    Qs <- diag(nrow(XX)) - tcrossprod(qr.Q(qr(XXs)))
    
    for (i in N_hat) {
      
      xx <- crossprod(sqrtW, XX[,i])
      beta_PIPE[i] <- t(xx) %*% (yy_pse - XXs %*% beta[c(TRUE, S)]) / crossprod(xx)
      sigma_PIPE[i] <- sqrt(1 / (t(xx) %*% Qs %*% xx))
      
    }
    
  }
  
  t <- beta_PIPE / sigma_PIPE
  pvalue <- (1 - pnorm(abs(t)))*2
  qvalue <- p.adjust(pvalue, method = "BH")
  
  if (fit$family != "gaussian") beta <- beta[-1]
  if (missing(sigma)) {sigma <- NA}
  
  if (posterior) {
    
    if (family == "gaussian") {
      ci <- ci_full_cond(beta_PIPE, lambda = lambda, se = sigma_PIPE,
                         alpha = 1 - level, penalty = penalty)
    } else {
      ci <- ci_full_cond(beta_PIPE, lambda = lambda, se = sigma_PIPE,
                         alpha = 1 - level, penalty = penalty, family = family, weights = weights)
    }
    
    
    lower <- ci$lower
    upper <- ci$upper
    
  } else {
    
    ci_width <- qnorm(1 - ((1 - level)/2)) * sigma_PIPE
    lower <- (beta_PIPE - ci_width)
    upper <- (beta_PIPE + ci_width)
    
  }
  
  res <- data.frame(
    variable = colnames(X),
    coef = beta / rescale_factorX,
    estimate = beta_PIPE / rescale_factorX,
    SE = sigma_PIPE,
    t = t,
    lower = lower / rescale_factorX,
    upper = upper / rescale_factorX,
    p.value = pvalue,
    p.adjust = qvalue,
    sigma = sigma,
    lamda = lambda,
    rescale_factor = rescale_factorX
  )
  
  return(res)
  
}
ci_full_cond <- function(z, lambda, se, alpha, penalty = "lasso", family = "gaussian", weights = NULL) {
  
  ## ADJUST LAMBDA FOR non gaussian
  if (family != "gaussian") {
    lambda_orig <- lambda
    lambda <- lambda / weights
  }
  
  ## Tails being transferred on to (log probability in each tail)
  obs_lw <- pnorm(0, z + lambda, se, log.p = TRUE)
  obs_up <- pnorm(0, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  
  obs_p_lw <- obs_lw + (z*lambda / se^2)
  obs_p_up <- obs_up - (z*lambda / se^2)
  
  ## Find the proportion of each to the overall probability
  frac_lw_log <- ifelse(is.infinite(exp(obs_p_lw - obs_p_up)), 0, obs_p_lw - obs_p_up - log(1 + exp(obs_p_lw - obs_p_up)))
  frac_up_log <- ifelse(is.infinite(exp(obs_p_up - obs_p_lw)), 0, obs_p_up - obs_p_lw - log(1 + exp(obs_p_up - obs_p_lw)))
  
  ps <- alpha / 2
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  lowers <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(log_one_minus_ps + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  )
  
  
  ps <- 1 - alpha / 2
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  uppers <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(log_one_minus_ps + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  )
  
  if (penalty %in% c("MCP", "SCAD")) {
    if (family != "gaussian") {lambda <- lambda_orig}
    uppers <- threshold_c(uppers, lambda, 3, family = family, weights = weights, penalty = penalty)
    lowers <- threshold_c(lowers, lambda, 3, family = family, weights = weights, penalty = penalty)
  }
  return(data.frame(lower = lowers, upper = uppers, lambda = lambda))
  
}
threshold_c <- function(bound, lambda, gamma, family = "gaussian", weights = NULL, penalty = "MCP") {
  
  # If the original code expected to never have bound == 0, we must ensure that condition.
  if (any(bound == 0)) {
    stop("An unexpected error occurred when adjusting the lasso penalty to MCP/SCAD.")
  }
  
  # Compute z_j depending on family
  if (family == "gaussian") {
    z_j <- bound
  } else {
    # Non-gaussian uses adaptive rescaling
    z_j <- bound * weights
  }
  
  # Add penalty adjustment
  z_j <- z_j + sign(z_j)*lambda
  
  if (penalty == "MCP") {
    # MCP thresholding
    # Condition for MCP
    cond <- abs(z_j) <= (gamma * lambda)
    
    # Apply MCP thresholding where condition holds
    z_j[cond] <- (gamma / (gamma - 1)) * soft_threshold(z_j[cond], lambda)
    
    # If non-gaussian, scale back by weights
    if (family != "gaussian") {
      z_j <- z_j / weights
    }
    return(z_j)
    
  } else if (penalty == "SCAD") {
    # SCAD thresholding
    cond1 <- abs(z_j) <= 2 * lambda
    cond2 <- (abs(z_j) > 2 * lambda) & (abs(z_j) <= gamma * lambda)
    
    out <- z_j
    # Region 1
    out[cond1] <- soft_threshold(z_j[cond1], lambda)
    
    # Region 2
    lambda_alt <- (gamma * lambda) / (gamma - 1)
    out[cond2] <- ((gamma - 1) / (gamma - 2)) * soft_threshold(z_j[cond2], lambda_alt)
    
    # Region 3: No change
    
    # If non-gaussian, scale back by weights
    if (family != "gaussian") {
      out <- out / weights
    }
    return(out)
    
  } else {
    stop("Unsupported penalty. Use 'MCP' or 'SCAD'.")
  }
}
soft_threshold <- function(z_j, lambda) {
  # Initialize an output vector of the same size as z_j
  out <- numeric(length(z_j))
  
  # Condition where z_j > lambda
  idx_pos <- z_j > lambda
  out[idx_pos] <- z_j[idx_pos] - lambda
  
  # Condition where |z_j| ≤ lambda
  idx_mid <- abs(z_j) <= lambda
  out[idx_mid] <- 0
  
  # Condition where z_j < -lambda
  idx_neg <- z_j < -lambda
  out[idx_neg] <- z_j[idx_neg] + lambda
  
  return(out)
}

