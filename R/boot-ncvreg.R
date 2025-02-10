#' Hybrid Bootstrap Confidence Intervals
#' 
#' Performs a hybrid bootstrapping approach to construct quantile based
#' confidence intervals around the original lasso/MCP/SCAD estimator.
#' Specifically, a traditional pairs bootstrap is performed with 1 adjustment:
#' if the bootstrap sample for a given covariate is zero, a random sample from
#' the full conditional posterior is used as the bootstrap sample instead.
#' This avoids the creation of intervals with endpoints exactly equal to zero.
#' 
#' The resulting intervals WILL NOT have exact nominal coverage for all
#' covariates. They are instead constructed in a way that overall coverage will
#' be approximately equal to nominal so long as the true distribution of betas
#' is Laplace and the covariates are independent. That said, in practice,
#' average coverage is fairly robust to these assumptions.
#' 
#' Note: Draws from the full conditional posterior are approximations for
#' MCP/SCAD or when \code{alpha} is not 1.
#' 
#' @param X       The design matrix, without an intercept. \code{boot_ncvreg}
#'                standardizes the data and includes an intercept by default.
#' @param y       The response vector.
#' @param fit     (optional) An object of class \code{ncvreg} or 
#'                \code{cv.ncvreg}. An object of class \code{ncvreg} simply
#'                provides data and penalty choices to \code{boot_ncvreg}. An
#'                object of class \code{cv.ncvreg} can in addition can provide
#'                information for selecting \code{lambda} and estimating
#'                \code{sigma2}. If provided, \code{y} should not be provided
#'                and \code{X} should only be provided if \code{fit} does not
#'                contain \code{X}.
#' @param lambda  (optional) The value of lambda to provide interval estimates
#'                for. If left missing will be selected using CV. If user wants
#'                to set the lambda sequence used to select \code{lambda} via
#'                cross validation, they should call \code{cv.ncvreg} separately
#'                and pass the resulting object to \code{fit}.
#' @param sigma2  (optional) The variance to use for the Hybrid sampling. If
#'                left missing will be set using the estimator suggested by Reid
#'                et. al. (2016) using CV.
#' @param cluster Bootstrapping and \code{cv.ncvreg} (if applicable) can be run
#'                in parallel across a cluster using the **parallel** package.
#'                The cluster must be set up in advance using the
#'                [parallel::makeCluster()] function from that package. The
#'                cluster must then be passed to \code{boot_ncvreg}.
#' @param seed    You may set the seed of the random number generator in order
#'                to obtain reproducible results. This is set for the overall
#'                process. If the user wishes to set a seed specifically for
#'                \code{cv.ncvreg} they should call it separately then pass the 
#'                fitted object as an argument to \code{fit}.
#' @param nboot   The number of bootstrap replications to use.
#' @param penalty The penalty to be applied to the model.  Either "lasso" (the
#'                default), "MCP", or "SCAD".
#' @param level   The confidence level required.
#' @param gamma   The tuning parameter of the MCP/SCAD penalty
#'                (see \code{ncvreg} for details). Default is 3 for MCP and 3.7
#'                for SCAD.
#' @param alpha   Tuning parameter for the Elastc net estimator which controls
#'                the relative contributions from the lasso/MCP/SCAD penalty and
#'                the ridge, or L2 penalty. `alpha=1` is equivalent to
#'                lasso/MCP/SCAD penalty, while `alpha=0` would be equivalent to
#'                ridge regression. However, `alpha=0` is not supported; `alpha`
#'                may be arbitrarily small, but not exactly 0.
#' @param returnCV    If \code{TRUE}, the \code{cv.ncvreg} fit will be returned
#'                    (if applicable).
#' @param return_boot If \code{TRUE}, the bootstrap draws will be returned.
#' @param verbose     If \code{TRUE}, messages are displayed indicating how
#'                    \code{lambda} and \code{sigma} are being selected.
#' @param ...         named arguments to be passed to \code{ncvreg} and
#'                    \code{cv.ncvreg}.
#'
#' @returns A list with:
#' \describe{
#'  \item{confidence_intervals}{A \code{data.frame} with the original point estimates along with lower and upper bounds of Hybrid CIs.}
#'  \item{lambda}{The value of \code{lambda} the \code{confidence_intervals} were constructed at.}
#'  \item{sigma2}{The value of \code{sigma2} used for the Hybrid bootstrap sampling.}
#'  \item{penalty}{The penalty the intervals correspond to.}
#'  \item{alpha}{The tuning parameter for the Enet estimator used.}
#' }
#' If a penalty other than "lasso" is used,
#' \describe{
#'  \item{gamma}{The tuning parameter for MCP/SCAD penalty.}
#' }
#' If \code{returnCV} is \code{TRUE} and a \code{cv.ncvreg} object was fit or supplied
#' \describe{
#'  \item{cv.ncvreg}{The \code{cv.ncvreg} fit used to estimate \code{lambda} and \code{sigma2} (if applicable).}
#' }
#' If \code{return_boot} is \code{TRUE}
#' \describe{
#'  \item{boot_draws}{A \code{data.frame} of the Hybrid bootstrap draws are returned.}
#' }
#' @export boot_ncvreg
#'
#' @examples
#' data(Prostate)
#' X <- Prostate$X
#' y <- Prostate$y
#' boot_ncvreg(X, y, level = 0.8)
boot_ncvreg <- function(X, y, fit, lambda, sigma2, cluster, seed,  nboot = 1000,
                        penalty = "lasso", level = 0.95,
                        gamma = switch(penalty, SCAD = 3.7, 3), alpha = 1,
                        returnCV = FALSE, return_boot = FALSE, verbose = FALSE,
                        ...) {
  
  if ((missing(X) | missing(y)) & (missing(fit) || !inherits(fit, c("cv.ncvreg", "ncvreg")))) {
    stop("Either X and y or a fit of class ncvreg or cv.ncvreg must be supplied.")
  }
  
  if (!missing(X) & missing(fit)) { # After &, just makes sure warning isn't duplicated as this is checked below as well
    if (any(attr(ncvreg::std(X), "scale") == 0)) warning("Some columns in X are singular. Intervals cannot be produced for corresponding covariates.")
  }
  
  cv_fit <- NULL
  if (!missing(fit)) {
    
    original_object_class <- class(fit)[1]
    
    if (inherits(fit, "cv.ncvreg")) {
      
      cv_fit <- fit
      fit <- cv_fit$fit
      
    }
    
    if (is.null(fit$X) & missing(X)) {
      stop(paste0("fit object missing X, please rerun ", original_object_class, " with returnX = TRUE or supply X directly along with the ", original_object_class, " object."))
    } 
    if (!missing(X) && is.null(attr(X, "scale"))) X <- std(X)
    if (!is.null(fit$X) && !missing(X) && !identical(fit$X, X)) {
      stop("X supplied along with ",  original_object_class, " object which also contains X and they are not the same. It is unclear which should be used.")
    } 
    if (missing(X)) X <- fit$X
    if (any(attr(X, "scale") == 0)) warning("Some columns in X are singular. Intervals cannot be produced for corresponding covariates.")
   
    if (!missing(y) && !identical(fit$y, y)) {
      stop(paste0("y supplied along with ",  original_object_class, " object which also contains y and they are not the same. It is unclear which should be used."))
    }
    if (missing(y)) y <- fit$y
    
    if (fit$family != "gaussian") {
      stop(paste0("fit of object class ", original_object_class, " was fit with ", fit$family, " family, but only 'gaussian' family is currently supported for boot_ncvreg"))
    }
    
  }
  
  # Coercion
  if (!is.matrix(X)) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  ## Seed handling
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  
  args <- list(...)
  cv.args <- args[names(args) %in% c("nfolds", "fold", "returnY", "trace")] ## returnY is depricated
  ncvreg.args <- args[names(args) %in% c("lambda.min", "nlambda", "eps", "max.iter", "dfmax", "warn")]
  
  if (length(cv.args) > 0 & !is.null(cv_fit)) {
    warning("Additional arguments for cv.ncvreg are ignored when a cv.ncvreg object is supplied")
  }
  if ("penalty.factor" %in% names(args)) {
    stop("Sorry, specification of alternative penality factors is not supported")
  }
  if (any(c("returnX", "convex") %in% names(args))) {
    message(paste0("Ignoring argument(s) ", paste0(names(args)[names(args) %in% c("returnX", "convex")], collapse = ", "), " they are set to FALSE in cv.ncvreg and any ncvreg objects fit are not accessible to user."))
  }
  
  if ("family" %in% names(args)) {
    stop("Ignoring argument 'family', only guassian family is currently supported")
  }
  
  original_coefs <- NULL
  if (is.null(cv_fit) & (missing(lambda) | missing(sigma2))) {
      ## Not needed unless we don't have sigma2 and lambda
      cv.args$X       <- X
      cv.args$y       <- y
      cv.args$penalty <- penalty
      cv.args$alpha   <- alpha
      cv.args$gamma   <- gamma
      if (!missing(cluster)) cv.args$cluster <- cluster
      cv_fit <- do.call("cv.ncvreg", c(cv.args, ncvreg.args))
      if (missing(fit)) fit <- cv_fit$fit ## Considering here if should just always use fit from cv_fit

  }
    
  if (verbose) {
    if (missing(lambda)) {
      message("Using cross validation to select lambda.")
    } else {
      message("Using user specified value for lambda.")
    }
    if (missing(sigma2)) {
      message("Estimating variance using Reid estimator.")
    } else {
      message("Using user specified value for sigma2.")
    }
  }
  
  if (missing(lambda)) lambda <- cv_fit$lambda.min
  if (missing(sigma2)) {
    ## Using fit from cv_fit here to ensure consistency
    yhat       <- cv_fit$fit$linear.predictors[,cv_fit$min]
    reid_coefs <- coef(cv_fit$fit, lambda = cv_fit$lambda.min)[-1]
    sh_lh      <- sum(reid_coefs != 0)
    sigma2     <- (max(length(y) - sh_lh, 1))^(-1) * sum((y - yhat)^2)
    
  }
  
  ## If X, y, lambda, sigma2 specified, still need an ncvreg fit
  if (missing(fit)) {
    ncvreg.args$X       <- X
    ncvreg.args$y       <- y
    ncvreg.args$penalty <- penalty
    ncvreg.args$alpha   <- alpha
    ncvreg.args$gamma   <- gamma
    fit <- do.call("ncvreg", ncvreg.args)
  }
  
  original_coefs <- coef(fit, lambda = lambda)[-1]
  
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "lambda", "sigma2", "ncvreg.args", "penalty", "gamma", "alpha"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    results <- parallel::parLapply(
      cl=cluster, X=1:nboot, fun=bootf, XX = X, yy=y, lambda = lambda,
      sigma2 = sigma2, ncvreg.args = ncvreg.args,
      penalty = penalty, alpha = alpha, gamma = gamma
    )
  }
  
  if (is.null(attr(X, "nonsingular"))) {
    ns <- 1:length(original_coefs)
  } else {
    ns <- attr(X, "nonsingular")
  }
  
  boot_draws <- matrix(nrow = nboot, ncol = length(original_coefs))
  for (i in 1:nboot) {
    if (!missing(cluster)) {
      boot_draws[i,ns] <- results[[i]]
    } else {
      boot_draws[i,ns] <- bootf(XX=X, yy=y, lambda = lambda, sigma2 = sigma2,
                   ncvreg.args = ncvreg.args,
                   penalty = penalty, alpha = alpha, gamma = gamma)
    }
  }
  
  colnames(boot_draws) <- names(original_coefs)
  cis <- compute_intervals(boot_draws, alpha = 1 - level)
  cis <- data.frame(
    estimates = original_coefs,
    cis
  )
  
  val <- list(
    confidence_intervals = cis,
    lambda               = lambda,
    sigma2               = sigma2,
    penalty              = penalty,
    alpha                = alpha
  )
  
  if (penalty != "lasso") val$gamma <- gamma
  if (return_boot) val$boot_draws <- boot_draws
  if (returnCV)    val$cv.ncvreg  <- cv_fit
  
  return(val)
  
}
bootf <- function(XX, yy, lambda, sigma2, ncvreg.args,
                  penalty = c("lasso", "MCP", "SCAD"),
                  alpha = 1, gamma = switch(penalty, SCAD = 3.7, 3),
                  ...) {
  
  p <- ncol(XX)
  n <- length(yy)
  
  boot_draws <- numeric(p)
  
  idx_new <- sample(1:n, replace = TRUE)
  ynew <- yy[idx_new]
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX[idx_new,,drop=FALSE])
  
  nonsingular <- attr(xnew, "nonsingular")
  p_nonsingular <- length(nonsingular)
  
  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale"))) {
    XXnonsingular <- attr(XX, "nonsingular")
    rescaleXX <- (attr(XX, "scale")[XXnonsingular][nonsingular])^(-1)
  } else {
    rescaleXX <- 1
  }
  full_rescale_factor <- rescale * rescaleXX
  
  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- setupLambda(
    xnew, ynew, "gaussian", alpha = alpha,
    nlambda, lambda.min = ifelse(n > p, .001, .05),
    penalty.factor=rep(1, ncol(xnew))
  )
  
  if (lambda < min(lambda_seq)) {
    
    message("Lambda too small, extending lambda sequence for bootstrap sample.")
    lambda_min <- lambda - (lambda / 100)
    lambda_seq <- setupLambda(
      xnew, ynew, "gaussian", alpha = alpha,
      lambda.min = lambda_min / max(lambda_seq), nlambda,
      penalty.factor=rep(1, ncol(xnew))
    )
    
  }
  
  if (lambda >= max(lambda_seq)) {
    message("Lambda too large, setting to lambda max for bootstrap sample.")
    lambda <- max(lambda_seq) - (max(lambda_seq) / 100)
  }
  
  ncvreg.args$X       <- xnew
  ncvreg.args$y       <- ynew
  ncvreg.args$lambda  <- lambda_seq
  ncvreg.args$penalty <- penalty
  ncvreg.args$alpha   <- alpha
  ncvreg.args$gamma   <- gamma
  
  ## Ignore user specified lambda.min, nlambda here since lambda sequence is being specified
  fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])
  modes <- coef(fit, lambda = lambda)[-1]
  
  ## For "elastic net" idea
  if (alpha < 1) {
    ynew <- c(ynew, rep(0, p_nonsingular))
    xnew <- rbind(xnew, sqrt(n*(1 - alpha)*lambda)*diag(p_nonsingular))
    xnew <- ncvreg::std(xnew)
    xnew <- xnew * sqrt(n / (n + p_nonsingular))
    lambda <- lambda * alpha
  }
  
  partial_residuals <- ynew - (
    as.numeric(xnew %*% modes) - (xnew * matrix(modes, nrow = nrow(xnew), ncol = ncol(xnew), byrow=TRUE))
  )
  z <- (1/n)*colSums(xnew * partial_residuals)
  if (sum(modes == 0) > 0) draws <- draw_full_cond(z[modes == 0], lambda, sigma2, n) ## Only where modes are 0
  
  if (penalty == "MCP") {
    draws <- sapply(draws, firm_threshold_c, lambda, gamma)
  } else if (penalty == "SCAD") {
    draws <- sapply(draws, scad_threshold_c, lambda, gamma)
  }
  
  if (sum(modes == 0) > 0) modes[modes == 0] <- draws
  boot_draws                          <- numeric(p)
  # boot_draws[nonsingular]             <- modes * full_rescale_factor
  tryCatch({
    boot_draws[nonsingular] <- modes * full_rescale_factor
  }, error = function(err) {
    message("Error occurred: ", err$message)
    message("Modes:")
    print(modes)
    message("Full rescale factor:")
    print(full_rescale_factor)
    message("Boot draws:")
    print(boot_draws)
    message("Nonsingular:")
    print(nonsingular)
  })
  
  boot_draws[!(1:p %in% nonsingular)] <- NA
  
  return(boot_draws)
  
}
draw_full_cond <- function(z, lambda, sigma2, n) {
  
  ## Tails being transferred on to (log probability in each tail)
  se <- sqrt(sigma2 / n)
  obs_lw <- pnorm(0, z + lambda, se, log.p = TRUE)
  obs_up <- pnorm(0, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  
  obs_p_lw <- obs_lw + ((z*lambda*n) / sigma2)
  obs_p_up <- obs_up - ((z*lambda*n) / sigma2)
  
  ## Find the proportion of each to the overall probability
  frac_lw_log <- ifelse(is.infinite(exp(obs_p_lw - obs_p_up)), 0, obs_p_lw - obs_p_up - log(1 + exp(obs_p_lw - obs_p_up)))
  frac_up_log <- ifelse(is.infinite(exp(obs_p_up - obs_p_lw)), 0, obs_p_up - obs_p_lw - log(1 + exp(obs_p_up - obs_p_lw)))
  
  ps <- runif(length(z))
  
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  draws <- ifelse(
    frac_lw_log >= log_ps,
    qnorm(log_ps + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
    qnorm(log_one_minus_ps + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
  )
  return(draws)
  
}
soft_threshold <- function(z_j, lambda) {
  
  if (z_j > lambda) {
    return(z_j - lambda)
  } else if (abs(z_j) <= lambda) {
    return(0)
  } else if (z_j < -lambda) {
    return(z_j + lambda)
  }
  
}
firm_threshold_c <- function(z_j, lambda, gamma) {
  
  z_j <- z_j + sign(z_j)*lambda
  
  if (abs(z_j) <= gamma*lambda) {
    return((gamma / (gamma - 1))*soft_threshold(z_j, lambda))
  } else {
    return(z_j)
  }
  
}
scad_threshold_c <- function(z_j, lambda, gamma) {
  
  z_j <- z_j + sign(z_j)*lambda
  
  if (abs(z_j) <= 2*lambda) {
    return(soft_threshold(z_j, lambda))
  } else if (abs(z_j) > 2*lambda & abs(z_j) <= gamma*lambda) {
    lambda_alt <- (gamma*lambda) / (gamma - 1)
    return(((gamma - 1) / (gamma - 2)) * soft_threshold(z_j, lambda_alt))
  } else {
    z_j
  }
  
}
compute_intervals <- function(draws, alpha = 0.2, quiet = FALSE) {
  
  any_nas <- any(as.logical(apply(draws, 2, function(x) sum(is.na(x)) > 0)))
  if (any_nas & !quiet) {
    warning("NAs in draws, this usually means a column in the bootstrapped X was singular.")
  }
  
  lowers <- apply(draws, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  uppers <- apply(draws, 2, function(x) quantile(x, 1 - (alpha / 2), na.rm = TRUE))
  
  return(data.frame(lower = lowers, upper = uppers))
  
}

