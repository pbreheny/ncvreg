#' Projection base test statistics and intervals
#' 
#' Constructs projection based test statistics that can be used to control FDR 
#' along with intervals for a penalized regression model.
#' 
#' The function constructs test statistics and intervals based off an
#' approximate projection onto the column space of the active features. The test
#' statistic can be used to control FDR and the intervals generally
#' have good coverage. However, both tend to be conservative with the
#' introduction of correlation.
#' 
#' The intervals produced can either be biased (like the point estimates) or 
#' debiased by setting the parameter \code{posterior} accordingly. The resulting
#' behavior is quite different. See references for more details.
#'
#' @param fit    An optional fit of class `ncvreg` or `cv.ncvreg`. If supplied,
#'               `X` should only be supplied if `fit` does not contain it.
#' @param lambda The penalty at which the tests and intervals are to
#'               be constructed. If left unspecified, will be selected using
#'               cross validation.
#' @param sigma Standard deviation estimate used to compute test statistic and 
#'              intervals. If left unspecified (default) it will be 
#'              estimated using the recommendation from Reid et al. (2016)
#' @param level the confidence level required.
#' @param posterior whether the intervals returned should be posterior intervals
#'              (default) or debiased intervals (if `FALSE`). Posterior
#'              intervals are constructed from distributions where the
#'              coefficient estimates are the is the posterior mode. Debiased
#'              intervals are constructed around the estimates.
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
#'   \item{estimate}{The debiased estimates.}
#'   \item{SE}{The standard errors.}
#'   \item{t}{The PIPE / Relaxed Lasso / LQA test statistics}
#'   \item{lower}{Interval lower bounds}
#'   \item{upper}{Intervals upper bounds}
#'   \item{p.value}{The unadjusted p-value}
#'   \item{p.adjust}{The Benhamini and Hochberg corrected p-value}
#'   \item{penalty}{The penalty used.}
#'   \item{lambda}{The lambda value the test statistics and intervals were constructed at.}
#'   \item{gamma}{The gamma value the test statistics and intervals were constructed at (for MCP/SCAD).}
#'   \item{alpha}{The alpha value the test statistics and intervals were constructed at.}
#'   \item{level}{The confidence level set for interval construction.}
#'   \item{sigma}{The standard deviation used for constructing the test statistis and intervals.}
#' }
#' 
#' @author Logan Harris, Patrick Breheny, and Biyue Dai
#'
#' @references
#' Harris L and Breheny P. (2025) A new perspective on high dimensional confidence intervals.
#' *arXiv preprint*, arXiv:2508.03504.
#' \url{https://arxiv.org/abs/2508.03504}
#' 
#' Dai B. (2019) Projection-based inference and model selection for penalized regression.
#' PhD dissertation, University of Iowa, Iowa City, IA.
#' \doi{10.17077/etd.005250}
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
    adj_num <- penalty_derivative(beta, lambda, alpha, penalty, gamma) 
    adjw <- diag(adj_num / abs(beta))
    
    ## Compute sigma if not provided
    if (missing(sigma)) {
      sigma <- sqrt(crossprod(yy - XX %*% beta) / (n - sum(beta != 0)))
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
    
    weights <- NULL
    
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
    adj_num <- penalty_derivative(beta[-1], lambda, alpha, penalty, gamma, weights)
    adjw <- diag(c(0, adj_num / (abs(beta[-1]) + 1e-9)))
    
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
    
    if (penalty == "lasso") {
        ci <- laplace_ci(z = beta_PIPE, lambda = lambda, se = sigma_PIPE,
                       alpha = 1 - level, gamma = gamma, enet_alpha = alpha,
                       weights = weights)
    } else {
        ci <- nonconvex_ci(z = beta_PIPE, lambda = lambda, se = sigma_PIPE,
                         alpha = 1 - level, gamma = gamma, penalty = penalty,
                         enet_alpha = alpha, weights = weights)
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
    SE = sigma_PIPE / rescale_factorX,
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
laplace_ci <- function(z, lambda, se, alpha, gamma, enet_alpha = 1, weights = NULL) {
  
  if (enet_alpha < 1) {
    
    z      <- z / (1 + (1-enet_alpha)*lambda)
    se     <- se / sqrt((1 + (1-enet_alpha)*lambda))
    lambda <- lambda * enet_alpha
    
  }
  
  if (!is.null(weights)) {
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
  
  return(data.frame(lower = lowers, upper = uppers, lambda = lambda))
  
}
## For LQA adjustment
penalty_derivative <- function(beta, lambda, alpha = 0.5,
                               penalty = c("lasso", "MCP", "SCAD"),
                               gamma = 3, weights = NULL) {
  ## work on the positive axis
  b <- abs(beta)
  penalty <- match.arg(penalty)
  
  ## ridge part
  d_ridge <- (1 - alpha) * b * lambda
  
  ## sparse‐penalty part
  # precompute the MCP and SCAD formulas
  if (!is.null(weights)) lambda  <- lambda / weights
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
logsumexp_vec <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}
log_norm_masses <- function(mu, s2, lambda, gamma) {
  
  # Derived quantities
  v_in   <- (gamma / (gamma - 1)) * s2
  sd_in  <- sqrt(v_in)
  v_out  <- s2
  sd_out <- sqrt(v_out)
  m_neg  <- (gamma / (gamma - 1)) * (mu + lambda)
  m_pos  <- (gamma / (gamma - 1)) * (mu - lambda)
  thr    <- gamma * lambda
  
  # Negative outside: N(mu, s2) over (-Inf, -thr)
  log_mass_neg_out <- pnorm(-thr, mean = mu, sd = sd_out, log.p = TRUE)
  
  # Negative inside: N(m_neg, v_in) over [-thr, 0)
  log_mass_neg_in <- log(
    pnorm(0, mean = m_neg, sd = sd_in) -
      pnorm(-thr, mean = m_neg, sd = sd_in)
  )
  
  # Positive inside: N(m_pos, v_in) over [0, thr]
  log_mass_pos_in <- log(
    pnorm(thr, mean = m_pos, sd = sd_in) -
      pnorm(0, mean = m_pos, sd = sd_in)
  )
  
  # Positive outside: N(mu, s2) over (thr, Inf)
  log_mass_pos_out <- pnorm(thr, mean = mu, sd = sd_out,
                            lower.tail = FALSE, log.p = TRUE)
  
  # Return in number line order
  c(
    neg_outside = log_mass_neg_out,
    neg_inside  = log_mass_neg_in,
    pos_inside  = log_mass_pos_in,
    pos_outside = log_mass_pos_out
  )
}
log_norm_masses_scad <- function(mu, s2, lambda, gamma) {

  # Cutpoints
  thr1 <- lambda
  thr2 <- gamma * lambda

  # Inner
  v_i  <- s2; sd_i <- sqrt(v_i)
  m_i_neg <- mu + lambda
  m_i_pos <- mu - lambda
  
  # Middle
  A   <- (gamma - 2) / (gamma - 1)
  v_m <- s2 / A; sd_m <- sqrt(v_m)
  cshift <- gamma / (gamma - 1)
  m_m_neg <- (mu + cshift * lambda) / A
  m_m_pos <- (mu - cshift * lambda) / A
  
  # Outside
  sd_o <- sqrt(s2)
  
  # Intervals
  log_neg_out <- pnorm(-thr2, mean = mu, sd = sd_o, log.p = TRUE)
  log_pos_out <- pnorm( thr2, mean = mu, sd = sd_o, lower.tail = FALSE, log.p = TRUE)
  
  log_neg_mid <- log(
    pnorm(-thr1, mean = m_m_neg, sd = sd_m) -
      pnorm(-thr2, mean = m_m_neg, sd = sd_m)
  )
  log_pos_mid <- log(
    pnorm( thr2, mean = m_m_pos, sd = sd_m) -
      pnorm( thr1, mean = m_m_pos, sd = sd_m)
  )
  
  log_neg_in <- log(
    pnorm(0,      mean = m_i_neg, sd = sd_i) -
      pnorm(-thr1,  mean = m_i_neg, sd = sd_i)
  )
  log_pos_in <- log(
    pnorm( thr1,  mean = m_i_pos, sd = sd_i) -
      pnorm(0,      mean = m_i_pos, sd = sd_i)
  )
  
  c(
    neg_outside = log_neg_out,
    neg_middle  = log_neg_mid,
    neg_inside  = log_neg_in,
    pos_inside  = log_pos_in,
    pos_middle  = log_pos_mid,
    pos_outside = log_pos_out
  )
}
log_constants <- function(mu, s2, lambda, gamma) {
  
  a <- (gamma - 1) / gamma  
  
  # Normalization constants from dnorm that we need to cancel
  logZ_out <- -0.5 * log(2 * pi * s2)
  logZ_in  <- -0.5 * log(2 * pi * (s2 / a))
  
  # Inside constants
  logC_neg <- (mu + lambda)^2 / (2 * a * s2) - mu^2 / (2 * s2) - logZ_in
  logC_pos <- (mu - lambda)^2 / (2 * a * s2) - mu^2 / (2 * s2) - logZ_in
  
  # Outside constant
  logC_out <- - (gamma * lambda^2) / (2 * s2) - logZ_out
  
  c(
    neg_outside = logC_out,
    neg_inside  = logC_neg,
    pos_inside  = logC_pos,
    pos_outside = logC_out
  )
}
log_constants_scad <- function(mu, s2, lambda, gamma) {
  
  # Cutpoints
  thr1 <- lambda
  thr2 <- gamma * lambda
  
  # Underlying normal variances per zone
  v_o <- s2                               
  v_i <- s2                               
  A   <- (gamma - 2) / (gamma - 1)        
  v_m <- s2 / A             
  
  # Remove dnorm normalization constants inside C
  logZ_o <- -0.5 * log(2 * pi * v_o)
  logZ_i <- -0.5 * log(2 * pi * v_i)
  logZ_m <- -0.5 * log(2 * pi * v_m)
  
  # Inner
  logC_ni <- ((mu + lambda)^2 - mu^2) / (2 * s2) - logZ_i   # neg_inside
  logC_pi <- ((mu - lambda)^2 - mu^2) / (2 * s2) - logZ_i   # pos_inside
  
  # Middle
  cshift   <- gamma / (gamma - 1)
  m_m_neg  <- (mu + cshift * lambda) / A
  m_m_pos  <- (mu - cshift * lambda) / A
  logC_nm  <- ( m_m_neg^2 / (2 * v_m) ) - (mu^2 / (2 * s2)) + (lambda^2 / (2 * (gamma - 1) * s2)) - logZ_m
  logC_pm  <- ( m_m_pos^2 / (2 * v_m) ) - (mu^2 / (2 * s2)) + (lambda^2 / (2 * (gamma - 1) * s2)) - logZ_m
  
  # Outside
  logC_out <- - ((gamma + 1) * lambda^2) / (2 * s2) - logZ_o
  
  c(
    neg_outside = logC_out,
    neg_middle  = logC_nm,
    neg_inside  = logC_ni,
    pos_inside  = logC_pi,
    pos_middle  = logC_pm,
    pos_outside = logC_out
  )
}
posterior_component_proportions <- function(mu, s2, lambda, gamma) {

  ln_vals <- log_norm_masses(mu, s2, lambda, gamma) +
    log_constants(mu, s2, lambda, gamma)
  
  # normalize in log-space to get proportions
  props <- exp(ln_vals - logsumexp_vec(ln_vals))
  
  # name + return
  names(props) <- c("neg_outside", "neg_inside", "pos_inside", "pos_outside")
  props
}
posterior_component_proportions_scad <- function(mu, s2, lambda, gamma) {
  
  ln_vals <- log_norm_masses_scad(mu, s2, lambda, gamma) +
    log_constants_scad(mu, s2, lambda, gamma)
  
  props  <- exp(ln_vals - logsumexp_vec(ln_vals))
  
  names(props) <- c("neg_outside","neg_middle","neg_inside","pos_inside","pos_middle","pos_outside")
  props
  
}
posterior_quantile <- function(mu, s2, lambda,
                               penalty = c("MCP","SCAD"),
                               gamma = switch(penalty, SCAD = 3.7, 3),
                               p) {
  
  penalty <- match.arg(penalty)
  
  ## Vectorize
  mu     <- as.numeric(mu)
  s2     <- as.numeric(s2)
  lambda <- as.numeric(lambda)
  p      <- as.numeric(p)
  
  lens <- c(length(mu), length(s2), length(lambda), length(p))
  n <- max(lens)
  
  if (length(mu)     == 1) mu     <- rep(mu,     n)
  if (length(s2)     == 1) s2     <- rep(s2,     n)
  if (length(lambda) == 1) lambda <- rep(lambda, n)
  if (length(p)      == 1) p      <- rep(p,      n)
  
  worker <- function(mu_i, s2_i, lam_i, p_i) {
    
    if (penalty == "MCP") {
      
      w       <- posterior_component_proportions(mu_i, s2_i, lam_i, gamma)
      cumw    <- cumsum(w)
      weights <- exp(log_norm_masses(mu_i, s2_i, lam_i, gamma) - log(w))
      
      thr    <- gamma * lam_i
      v_in   <- (gamma/(gamma - 1)) * s2_i
      sd_in  <- sqrt(v_in)
      sd_out <- sqrt(s2_i)
      m_neg  <- (gamma/(gamma - 1)) * (mu_i + lam_i)
      m_pos  <- (gamma/(gamma - 1)) * (mu_i - lam_i)
      
      k <- which(p_i <= cumw)[1]
      prev <- if (k == 1) 0 else cumw[k-1]
      
      if (k == 1) {
        q <- qnorm(p_i * weights[1], mean = mu_i, sd = sd_out)
      } else if (k == 2) {
        p_loc <- weights[2]*(p_i - prev) + pnorm(-thr, m_neg, sd_in)
        q <- qnorm(p_loc, mean = m_neg, sd = sd_in)
      } else if (k == 3) {
        p_loc <- weights[3]*(p_i - prev) + pnorm(0, m_pos, sd_in)
        q <- qnorm(p_loc, mean = m_pos, sd = sd_in)
      } else {
        p_loc <- weights[4]*(p_i - prev) + pnorm(thr, mu_i, sd_out)
        q <- qnorm(p_loc, mean = mu_i, sd = sd_out)
      }
      return(q)
    }
    
    ## SCAD
    w <- posterior_component_proportions_scad(mu_i, s2_i, lam_i, gamma)
    cumw    <- cumsum(w)
    weights <- exp(log_norm_masses_scad(mu_i, s2_i, lam_i, gamma) - log(w))
    
    thr1 <- lam_i
    thr2 <- gamma * lam_i
    
    # inner
    v_i   <- s2_i
    sd_i  <- sqrt(v_i)
    m_i_neg <- mu_i + lam_i
    m_i_pos <- mu_i - lam_i
    
    # middle
    A      <- (gamma - 2) / (gamma - 1)
    v_m    <- s2_i / A
    sd_m   <- sqrt(v_m)
    cshift <- gamma / (gamma - 1)
    m_m_neg <- (mu_i + cshift * lam_i) / A
    m_m_pos <- (mu_i - cshift * lam_i) / A
    
    # outside
    sd_o <- sqrt(s2_i)
    
    k <- which(p_i <= cumw)[1]
    prev <- if (k == 1) 0 else cumw[k-1]
    
    if (k == 1) {
      q <- qnorm(p_i * weights[1], mean = mu_i, sd = sd_o)
    } else if (k == 2) {
      p_loc <- weights[2]*(p_i - prev) + pnorm(-thr2, mean = m_m_neg, sd = sd_m)
      q <- qnorm(p_loc, mean = m_m_neg, sd = sd_m)
    } else if (k == 3) {
      p_loc <- weights[3]*(p_i - prev) + pnorm(-thr1, mean = m_i_neg, sd = sd_i)
      q <- qnorm(p_loc, mean = m_i_neg, sd = sd_i)
    } else if (k == 4) {
      p_loc <- weights[4]*(p_i - prev) + pnorm(0,       mean = m_i_pos, sd = sd_i)
      q <- qnorm(p_loc, mean = m_i_pos, sd = sd_i)
    } else if (k == 5) {
      p_loc <- weights[5]*(p_i - prev) + pnorm(thr1,    mean = m_m_pos, sd = sd_m)
      q <- qnorm(p_loc, mean = m_m_pos, sd = sd_m)
    } else {
      p_loc <- weights[6]*(p_i - prev) + pnorm(thr2,    mean = mu_i,    sd = sd_o)
      q <- qnorm(p_loc, mean = mu_i, sd = sd_o)
    }
    q
  }
  
  vapply(seq_len(n), function(i) worker(mu[i], s2[i], lambda[i], p[i]), numeric(1))
}
nonconvex_ci <- function(z, se, lambda, gamma, penalty, alpha, enet_alpha = 1,
                         weights = NULL) {
  
  if (enet_alpha < 1) {
    z      <- z / (1 + (1-enet_alpha)*lambda)
    se     <- se / sqrt((1 + (1-enet_alpha)*lambda))
    lambda <- lambda * enet_alpha
  }
  
  if (!is.null(weights)) {
    lambda <- lambda / weights
  }
  
  lowers <- posterior_quantile(z, se^2, lambda, penalty, gamma, alpha / 2)
  uppers <- posterior_quantile(z, se^2, lambda, penalty, gamma, 1 - alpha / 2)
  
  return(data.frame(lower = lowers, upper = uppers, lambda = lambda))
  
}