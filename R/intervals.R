#' Projection base test statistic and confidence intervals
#' 
#' Constructs projection based test statistics that can be used to control FDR 
#' along with confidence intervals for a penalized regression model.
#' 
#' The function constructs test statistics and confidence intervals based off an
#' approximate projection onto the column space of the active features. The test
#' statistic can be used to control FDR and the confidence intervals generally
#' have good coverage. However, both tend to be conservative with the
#' introduction of correlation.
#'
#' @param fit    An optional fit of class `ncvreg` or `cv.ncvreg`. If supplied,
#'               `X` should only be supplied if `fit` does not contain it.
#' @param lambda The penalty at which the tests and confidence intervals are to
#'               be constructed. If left unspecified, will be selected using
#'               cross validation.
#' @param sigma Standard deviation estimate used to compute test statistic and 
#'              confidence intervals. If left unspecified (default) it will be 
#'              estimated using the recommendation from Reid et al. (2016)
#' @param level the confidence level required.
#' @param posterior whether the intervals returned should be posterior intervals
#'              (default) or debiased intervals (if `FALSE`). Posterior
#'              intervals are constructed from distributions where the estimate
#'              is the posterior mode. Debiased intervals are constructed
#'              around the test statistic.
#' @param relaxed whether the relaxed lasso based statistic / intervals should
#'              be used. Default is `FALSE` in which case PIPE based intervals
#'              are constructed (recommended). This affects the estimate.
#' @param adjust_projection whether a Local Quadratic Approximation should be
#'              used in the projection for determining variance of the PIPE 
#'              based test statistic. Default is `FALSE` as more research has
#'              been done without this adjustment, however without this
#'              adjustment, statistics and intervals may be over conservative in
#'              the presence of correlation. This affects the SE.
#' @param X     The original design matrix supplied to `fit`. Required if `fit`
#'              does not contain `X`.
#'
#' @returns An `data.frame` containing the following columns:
#' \describe{
#'   \item{variable}{`colnames(X)`}
#'   \item{coef}{The original estimates at the specified parameters (`lambda`, `gamma`, `alpha`)}
#'   \item{estimate}{The PIPE / Relaxed Lasso estimates}
#'   \item{SE}{The PIPE / LQA standard errors. The Relaxed Lasso SEs are the same as PIPE's.}
#'   \item{t}{The PIPE / Relaxed Lasso / LQA test statistics}
#'   \item{lower}{CI lower bounds}
#'   \item{upper}{CI upper bounds}
#'   \item{p.value}{The unadjusted p-value}
#'   \item{p.adjust}{The Benhamini and Hochberg corrected p-value}
#'   \item{penalty}{The penalty used.}
#'   \item{lambda}{The lambda value the test statistics and CIs were constructed at.}
#'   \item{gamma}{The gamma value the test statistics and CIs were constructed at (for MCP/SCAD).}
#'   \item{alpha}{The alpha value the test statistics and CIs were constructed at.}
#'   \item{level}{The confidence level set for CI construction.}
#'   \item{sigma}{The standard deviation used for constructing the test statistis and CIs.}
#' }
#' 
#' @author Logan Harris, Patrick Breheny, and Biyue Dai
#' 
#' @examples
#' # Linear regression (SCAD-Net penalty, PIPE intervals, pass ncvreg object)
#' fit <- ncvreg(Prostate$X, Prostate$y, penalty = "SCAD", alpha = 0.9)
#' intervals(fit) |> head()
#' 
#' # Logistic regression (lasso penalty, LQA intervals, pass cv.ncvreg object) 
#' data(Heart)
#' cv_fit <- cv.ncvreg(Heart$X, Heart$y, family="binomial", penalty = "lasso")
#' intervals(cv_fit, adjust_projection = TRUE) |> head()
#' 
#' @export intervals
intervals <- function(fit, lambda, sigma, level = 0.95,
                      posterior = TRUE, relaxed = FALSE, 
                      adjust_projection = FALSE, X = NULL
) {
  
  if (!inherits(fit, c("cv.ncvreg", "ncvreg"))) {
    stop("fit must be of class ncvreg or cv.ncvreg")
  }
    
  original_object_class <- class(fit)[1]
  
  if (inherits(fit, "cv.ncvreg")) {
    
    cv_fit <- fit
    fit <- cv_fit$fit
    
    if (missing(lambda)) lambda <- cv_fit$lambda.min
    if (lambda == max(cv_fit$lambda)) {lambda <- lambda * 0.999}
    if (lambda == min(cv_fit$lambda)) {lambda <- lambda * 1.001}
    
  }
  
  if (is.null(fit$X) & is.null(X)) {
    stop(paste0(
      "fit object missing X, please rerun ",
      original_object_class,
      " with returnX = TRUE or supply X directly along with the ",
      original_object_class,
      " object."
    ))
  }
  if (!missing(X) && !is.null(fit$X)) {
    if (is.null(attr(X, "scale"))) {
      Xcheck <- std(X)
    } else {
      Xcheck <- X
    }
    if (!identical(fit$X, Xcheck)) {
      stop(
        "X supplied along with ",
        original_object_class,
        " object which also contains X and they are not the same.",
        "It is unclear which should be used."
      )
    } 
  }
  
  if (any(fit$penalty.factor != 1)) {
    stop("Alternate penalty factors are currently not supported.")
  }
  
  ## Use the values from ncvreg or cv.ncvreg object
  penalty <- fit$penalty
  gamma <- fit$gamma
  alpha <- fit$alpha
  family <- fit$family
  
  ## Get standardized X and its scale for later rescaling
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
  if (any(attr(XX, "scale") == 0)) {
    stop(
      "Some columns in X are singular.",
      "Intervals cannot be produced for corresponding covariates.",
      "Please remove these columns before preceeding."
    )
  }
  
  # Get y and center if gaussian
  yy <- fit$y
  if (fit$family == "gaussian") yy <- yy - mean(yy)
  
  p <- ncol(XX)
  n <- nrow(XX)
  
  ## Get lambda if not yet defined
  if (missing(lambda)) {
    ## Select value for lambda if not provided or not extracted from cv_fit
    cv_fit <- cv.ncvreg(
      XX, yy, family = family, penalty = penalty,
      gamma = gamma, alpha = alpha, returnX = TRUE
    ) 
    lambda <- cv_fit$lambda.min 
    if (lambda == max(cv_fit$lambda)) {lambda <- lambda * 0.999}
    if (lambda == min(cv_fit$lambda)) {lambda <- lambda * 1.001}
  } 
  
  # initialize
  beta_PIPE <- numeric(p)
  sigma_PIPE <- numeric(p)
  
  if (fit$family == "gaussian"){
    
    # Load model with standardized coefficients
    beta <- coef(fit, lambda = lambda)[-1]
    beta <- beta * rescale_factorX
    
    ## Get S_hat
    S <- beta != 0
    S_hat <- which(S)
    N_hat <- which(!S)
    
    ## For LQA adjustment
    ## alpha treated as 1 here because of how this fits with the augmentation
    adj_num <- penalty_derivative(beta, lambda, 1, penalty, gamma) 
    adjw <- diag(adj_num / abs(beta))
    
    ## Compute sigma if not provided
    if (missing(sigma)) {
      sigma <- sqrt(crossprod(yy - XX %*% beta) / (n - sum(beta != 0)))
    }
    
    ## NET data augmentation approach
    if (alpha < 1 & posterior) {
      yy <- c(yy, rep(0, p))
      XX <- rbind(XX, sqrt(n*(1 - alpha)*lambda)*diag(p))
      lambda <- lambda * alpha
    }
    
    ## Compute PIPE test statistic
    yhat <- as.numeric(XX %*% beta)
    partial_residuals <- (yy - yhat) + 
      (XX * matrix(beta, nrow = nrow(XX), ncol = ncol(XX), byrow = TRUE))
    beta_PIPE <- (1/n) * colSums(XX * partial_residuals)
    
    ## Compute PIPE Variance for features in the estimated support
    for (i in S_hat){
      
      S_hat_i <- S
      S_hat_i[i] <- FALSE
      
      ## Compute variance
      if (sum(S_hat_i) > 0) {
        Xsi <- XX[,S_hat_i,drop = FALSE]
        if (!adjust_projection) {
          Qsi <- diag(nrow(XX)) - tcrossprod(qr.Q(qr(Xsi)))
        } else {
          adji <- adjw[S_hat_i, S_hat_i]
          Qsi <- diag(nrow(XX)) - 
            Xsi %*% solve(((1/n)*crossprod(Xsi) + adji)) %*% ((1/n) * t(Xsi))
        }
        adjusted_n <- t(XX[,i,drop = FALSE]) %*% Qsi %*% XX[,i,drop = FALSE]
      } else {
        adjusted_n <- nrow(XX)
        Qsi <- diag(nrow(XX))
      }
      
      sigma_PIPE[i] <- sqrt(sigma^2 / adjusted_n)
      
      ## If relaxed, update estimates
      if (relaxed) {
        beta_PIPE[i] <- (t(XX[,i,drop=FALSE]) %*% Qsi %*% yy) / adjusted_n
      }
      
    }
    
    # Same for null features
    adjs <- adjw[S_hat,S_hat]
    if (length(S_hat) > 0 & length(N_hat) > 0) {
      Xs <- XX[,S_hat,drop=FALSE]
      if (!adjust_projection) {
        Qs <- diag(nrow(XX)) - tcrossprod(qr.Q(qr(Xs)))
      } else {
        Qs <- diag(nrow(XX)) -
          Xs %*% solve(((1/n)*crossprod(Xs) + adjs)) %*% ((1/n) * t(Xs))
      }
    } else {
      Qs <- diag(nrow(XX))
    }
    
    for (i in N_hat) {
      
      adjusted_n <- t(XX[,i,drop = FALSE]) %*% Qs %*% XX[,i,drop = FALSE]
      sigma_PIPE[i] <- sqrt(sigma^2 / adjusted_n)
      
      if (relaxed) {
        beta_PIPE[i] <- (t(XX[,i,drop=FALSE]) %*% Qs %*% yy) / adjusted_n
      }
      
    }
    
  } else {
    
    ## For poisson and logistic penalized regression
    ## Refit the model with standardized X, I don't love this... 
    ## but so far have not been able to find an alternate solution
    if (!all(abs(attr(XX, "scale") - 1) < 1e-8)) {
      fit <- ncvreg::ncvreg(
        XX, yy, family = fit$family, penalty = fit$penalty,
        gamma = fit$gamma, alpha = fit$alpha
      )
    }
    
    ## Get estimates and predicted values
    pii <- predict(fit, XX, type = "response", lambda = lambda)
    beta <- coef(fit, lambda = lambda)
    S <- beta[-1] != 0
    S_hat <- which(S)
    N_hat <- which(!S)
    
    ## Construct response specific components
    if (fit$family == "binomial") {
      A <- pii * (1 - pii)
    } else if (fit$family == "poisson") {
      A <- pii
    }
    
    v <- yy - pii
    W <- diag(A)
    sqrtW <- diag(A^(1/2))
    Y_pse <- diag(1/A) %*% v + predict(fit, XX, type = "link", lambda = lambda)
    yy <- sqrtW %*% Y_pse
    X_int <- sqrtW %*% cbind(rep(1, nrow(XX)), XX)
    XX <- sqrtW %*% XX
    weights <- colSums(XX^2) / n
    
    ## For LQA approach
    adj_num <- penalty_derivative(beta[-1], lambda, 1, penalty, gamma)
    adjw <- diag(c(0, adj_num / (abs(beta[-1]) + 1e-9)))
    
    ## For NET Augmentation
    if (alpha < 1 & posterior) {
      yy <- c(yy, rep(0, p))
      XX <- rbind(XX, sqrt(n*(1 - alpha)*lambda)*diag(p))
      X_int <- rbind(
        X_int,
        cbind(rep(0, p), sqrt(n*(1 - alpha)*lambda)*diag(p))
      )
      lambda <- lambda * alpha
    }
    
    # Compute pipe for features in the estimated support
    for (i in S_hat) {
      
      S_hat_i <- S
      S_hat_i[i] <- FALSE
      
      Xsi <- X_int[, c(TRUE, S_hat_i), drop = FALSE]
      if (sum(S_hat_i) > 0) {
        if (!adjust_projection) {
          Qsi <- diag(nrow(Xsi)) - tcrossprod(qr.Q(qr(Xsi)))
        } else {
          adji <- adjw[c(TRUE, S_hat_i), c(TRUE, S_hat_i)]
          Qsi <- diag(nrow(XX)) -
            Xsi %*% solve((1/n)*crossprod(Xsi) + adji) %*% ((1/n) * t(Xsi))
        }
      } else {
        Qsi <- diag(nrow(XX))
      }
      
      ## Compute the statistic
      if (relaxed) {
        beta_PIPE[i] <- (t(XX[,i,drop=FALSE]) %*% Qsi %*% yy) / 
          (t(XX[,i,drop=FALSE]) %*% Qsi %*% XX[,i,drop=FALSE])
      } else {
        beta_PIPE[i] <- 
          (t(XX[,i,drop=FALSE]) %*% (yy - Xsi %*% beta[c(TRUE, S_hat_i)])) /
          (t(XX[,i,drop=FALSE]) %*% XX[,i,drop=FALSE])
      }
      
      ## Compute SE
      sigma_PIPE[i] <- sqrt(
        1 / (t(XX[,i,drop=FALSE]) %*% Qsi %*% XX[,i,drop=FALSE])
      )
      
    }
    
    # Compute pipe for null features
    Xs <- X_int[,c(TRUE, S), drop = FALSE]
    adjs <- adjw[c(TRUE, S),c(TRUE, S)]
    if (length(S_hat) > 0 & length(N_hat) > 0) {
      if (!adjust_projection) {
        Qs <- diag(nrow(Xs)) - tcrossprod(qr.Q(qr(Xs)))
      } else {
        Qs <- diag(
          nrow(XX)) -
          Xs %*% solve(((1/n)*crossprod(Xs) + adjs)) %*% ((1/n) * t(Xs)
        )
      }
    } else {
      Qs <- diag(nrow(XX))
    }
    
    for (i in N_hat) {
      if (relaxed) {
        beta_PIPE[i] <- (t(XX[,i,drop=FALSE]) %*% Qs %*% yy) /
          (t(XX[,i,drop=FALSE]) %*% Qs %*% XX[,i,drop=FALSE])
      } else {
        beta_PIPE[i] <- 
          (t(XX[,i,drop=FALSE]) %*% (yy - Xs %*% beta[c(TRUE, S)])) /
          (t(XX[,i,drop=FALSE]) %*% XX[,i,drop=FALSE])
      }
      
      sigma_PIPE[i] <- sqrt(
        1 / (t(XX[,i,drop=FALSE]) %*% Qs %*% XX[,i,drop=FALSE])
      )
      
    }
    
  }
  
  ## Compute test statistic, pvalue and qvalue
  t <- beta_PIPE / sigma_PIPE
  pvalue <- (1 - pnorm(abs(t)))*2
  qvalue <- p.adjust(pvalue, method = "BH")
  
  ## If poisson / logistic remove intercept before returning
  if (fit$family != "gaussian") beta <- beta[-1]
  
  ## If missing sigma (poisson / logistic) set sigma to missing
  if (missing(sigma)) {sigma <- NA}
  
  ## Compute intervals
  if (posterior) {
    
    if (family == "gaussian") {
      
      ci <- ci_full_cond(z = beta_PIPE, lambda = lambda, se = sigma_PIPE,
                         alpha = 1 - level, gamma = gamma, penalty = penalty)
      
    } else {
      
      ci <- ci_full_cond(z = beta_PIPE, lambda = lambda, se = sigma_PIPE,
                         alpha = 1 - level, gamma = gamma, penalty = penalty,
                         family = family, weights = weights)
      
    }
    
    lower <- ci$lower
    upper <- ci$upper
    
  } else {
    
    ci_width <- qnorm(1 - ((1 - level)/2)) * sigma_PIPE
    lower <- (beta_PIPE - ci_width)
    upper <- (beta_PIPE + ci_width)
    
  }
  
  ## Return results
  res <- data.frame(
    variable = names(beta),
    coef = beta / rescale_factorX,
    estimate = beta_PIPE / rescale_factorX,
    SE = sigma_PIPE,
    t = t,
    lower = lower / rescale_factorX,
    upper = upper / rescale_factorX,
    p.value = pvalue,
    p.adjust = qvalue,
    penalty = penalty,
    lambda = lambda,
    gamma  = gamma,
    alpha  = alpha,
    level  = level,
    sigma = sigma
  )
  
  return(res)
  
}
## Function for computing posterior intervals
ci_full_cond <- function(z, lambda, se, alpha, gamma, penalty = "lasso",
                         family = "gaussian", weights = NULL) {
  
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
  frac_lw_log <- ifelse(
    is.infinite(exp(obs_p_lw - obs_p_up)),
    0,
    obs_p_lw - obs_p_up - log(1 + exp(obs_p_lw - obs_p_up))
  )
  frac_up_log <- ifelse(
    is.infinite(exp(obs_p_up - obs_p_lw)),
    0,
    obs_p_up - obs_p_lw - log(1 + exp(obs_p_up - obs_p_lw))
  )
  
  ps <- alpha / 2
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  lowers <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(
      log_one_minus_ps + obs_up - frac_up_log, z - lambda, se,
      lower.tail = FALSE, log.p = TRUE
    )
  )
  
  ps <- 1 - alpha / 2
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  uppers <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(
      log_one_minus_ps + obs_up - frac_up_log, z - lambda, se,
      lower.tail = FALSE, log.p = TRUE
    )
  )
  
  if (penalty %in% c("MCP", "SCAD")) {
    ## Return lambda to original for MCP/SCAD approximation
    if (family != "gaussian") {lambda <- lambda_orig} 
    ## Do adjustment for MCP / SCAD
    uppers <- threshold_c(
      uppers, lambda, gamma, family = family,
      weights = weights, penalty = penalty
    )
    lowers <- threshold_c(
      lowers, lambda, gamma, family = family,
      weights = weights, penalty = penalty
    )
  }
  
  return(data.frame(lower = lowers, upper = uppers, lambda = lambda))
  
}
threshold_c <- function(bound, lambda, gamma, family = "gaussian", 
                        weights = NULL, penalty = "MCP") {
  
  # Compute z_j depending on family
  if (family == "gaussian") {
    z_j <- bound
  } else {
    z_j <- bound * weights ## Return bounds to original scal
  }
  
  # Add penalty adjustment, i.e. invert soft threshold
  z_j <- z_j + sign(z_j) * lambda
  
  ## Re-apply appropriate penalty
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
    out[cond2] <- ((gamma - 1) / (gamma - 2)) *
      soft_threshold(z_j[cond2], lambda_alt)
    
    # Region 3: No change
    
    # If non-gaussian, scale back by weights
    if (family != "gaussian") {
      out <- out / weights
    }
    return(out)
    
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
## For LQA adjustment
penalty_derivative <- function(beta, lambda, alpha = 0.5,
                               penalty = c("lasso", "MCP", "SCAD"),
                               gamma = 3) {
  ## work on the positive axis
  b <- abs(beta)
  penalty <- match.arg(penalty)
  
  ## ridge part
  d_ridge <- (1 - alpha) * b * lambda
  
  ## sparse‐penalty part
  # precompute the MCP and SCAD formulas
  mu_mcp  <- pmax(lambda - b/gamma, 0)
  mu_scad <- pmin(lambda, pmax((gamma*lambda - b)/(gamma - 1), 0))
  
  # then pick by penalty in one nested ifelse
  d_sparse_raw <- ifelse(
    penalty == "lasso", lambda,
    ifelse(
      penalty == "MCP",  mu_mcp,
      ifelse(
        penalty == "SCAD", mu_scad,
        NA_real_
      )
    )
  )
  
  # combine ridge + sparse
  d_ridge + alpha * d_sparse_raw
}