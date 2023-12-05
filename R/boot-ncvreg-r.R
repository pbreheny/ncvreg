#' Bootstrap for ncvreg with lasso penalty
#' 
#' Perform bootstrapping for the lasso at a single value of the regularization
#' parameter lambda.
#' 
#' At a specified value of lambda and variance 
#' (which are selected / estimated by default if not provided)
#' the code does \code{nboot} iterations of taking a pairs bootstrap sample of
#' \code{X} and \code{y}, fitting \code{ncvreg}, and then uses this fit to
#' obtain \code{significance_level * 100}% central credible intervals from 
#' marginal posteriors based on a Laplace prior and a Normal likelihood 
#' specified with partial residuals.
#'
#' @param X The design matrix, without an intercept, as in \code{ncvreg}
#' @param y The response vector, as in \code{ncvreg}
#' @param cv_fit Instead of \code{X} and \code{y}, an object of type 
#' \code{cv.ncvreg} with \code{returnX = TRUE} and \code{penalty = "lasso"}
#' @param lambda A positive scalar value of type numeric, if left unspecified, selected using \code{cv.ncvreg}
#' @param sigma2 A known / estimated value for the variance used in the Laplace prior and Normal likelihood, if left unspecified, selected using \code{cv.ncvreg} only if \code{lambda} is also left as unspecified
#' @param significance_level Desired confidence level for determining interval width, between 0 and 1. Default is 0.8.
#' @param nboot The number of bootstrap iterations. Default is 100.
#' @param ... 
#' @param cluster \code{cv.ncvreg} and \code{cv.ncvsurv} can be run in parallel
#' across a cluster using the \code{parallel} package.  The cluster must be set
#' up in advance using the \code{makeCluster} function from that package.  The
#' cluster must then be passed to \code{cv.ncvreg}/\code{cv.ncvsurv} (see
#' example).
#' @param seed You may set the seed of the random number generator in order to
#' obtain reproducible results. Note that the seed is set at the start of the 
#' code an thus applies to \code{cv.ncvreg} as well as the bootstrap if \code{cv.ncvreg}
#' is run to select \code{lambda} and estimate \code{sigma2}. If the user would
#' like to specify a different seed to \code{cv.ncvreg}, they should manually
#' specify the call to \code{cv.ncvreg} in the arguments.
#' @param returnCV Whether to return \code{cv.ncvreg} object, \code{FALSE} by default.
#'
#' @return An object with S3 class \code{boot.ncvreg} containing: \describe { 
#' \item{lowers}{A \code{nboot} by \code{ncol(X)} matrix of the lower bounds
#' from each bootstrap iteration for each variable.} \item{uppers}{A 
#' \code{nboot} by \code{ncol(X)} matrix of the upper bounds
#' from each bootstrap iteration for each variable.} \item{modes}{A 
#' \code{nboot} by \code{ncol(X)} matrix of the posterior modes (lasso solution)
#' from each bootstrap iteration for each variable.} \item{cv.ncvreg}{If 
#' \code{returnCV == TRUE}, an object of type \code{cv.ncvreg} that was fit to
#' \code{X} and \code{y} to optionally select \code{lambda} and estimate 
#' \code{sigma2}.}}
#' @seealso \code{\link{ncvreg}}
#' 
#' @examples
#' 
#' library(hdrm)
#' dat <- readData(whoari)
#' 
#' tmp <- boot.ncvreg(dat$X, dat$y)
#' 
#' ## Specifying cv_fit
#' tmp <- boot.ncvreg(cv_fit = cv.ncvreg(dat$X, dat$y, penalty = "lasso", returnX = TRUE))
#' 
#' @export boot.ncvreg.r
boot.ncvreg.r <- function(X, y, cv_fit, lambda, sigma2, significance_level = 0.8, nboot = 100, ..., cluster, seed, returnCV=FALSE, verbose = TRUE, time = FALSE, quantiles = "sample") {
  
  if (time) tic(msg = "Overall")
  
  if (time) tic(msg = "Checks")
  if ((missing(X) | missing(y)) & (missing(cv_fit) || class(cv_fit) != "cv.ncvreg")) {
    stop("Either X and y or an object of class cv.ncvreg must be supplied.")
  }
  
  if (!missing(cv_fit)) {
    if (!all(c("X", "y") %in% names(cv_fit$fit))) {
      stop("fit object in cv_fit missing X and y, please rerun cv.ncvreg with returnX = TRUE or supply X and y directly (without also specifying a cv.ncvreg object)") 
    }
    
    if(cv_fit$fit$penalty != "lasso") {
      stop(paste0("cv.ncvreg fit with ", cv_fit$fit$penalty, " penalty, but only 'lasso' penalty is currently supported for boot.ncvreg"))
    }
    
    if(cv_fit$fit$family != "gaussian") {
      stop(paste0("cv.ncvreg fit with ", cv_fit$fit$family, " family, but only 'gaussian' family is currently supported for boot.ncvreg"))
    }
  }
  
  if (!missing(y) & !missing(cv_fit)) {
    warning("Ignoring supplied values of y, using y from supplied cv.ncvreg object")
  }
  
  if (!missing(X) & !missing(cv_fit)) {
    if ("X" %in% names(cv_fit$fit)) {
      warning("Ignoring supplied values of X, using X from supplied cv.ncvreg object") 
    }
  }
  
  if (!missing(X)) {
    rescale_original <- FALSE
  } else {
    rescale_original <- TRUE
  }
  
  # Coercion
  if (!missing(X)) {
    if (!is.matrix(X)) {
      tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (storage.mode(X)=="integer") storage.mode(X) <- "double" 
  }
  if (!missing(y)) {
    if (!is.double(y)) {
      tmp <- try(y <- as.double(y), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
    } 
  }
  
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  
  args <- list(...)
  if (length(args) > 0 & !missing(cv_fit)) {
    warning("Additional arguments are ignored when cv.ncvreg object supplied")
  }
  
  if (length(args) > 0 & any(is.null(names(args)))) {
    stop("Please supply names for all additional arguments passed to ...")
  }
  
  
  cv.args <- args[names(args) %in% c("nfolds", "fold", "returnY", "trace")]
  ncvreg.args <- args[names(args) %in% c("lambda.min", "nlambda", "eps", "max.iter", "dfmax")]
  if ("penalty.factor" %in% names(args)) {
    stop("Sorry, specification of alternative penality factors is not yet supported")
  }
  
  if (any(c("returnX", "warn", "convex") %in% names(args))) {
    warning(paste0("Ignoring argument(s) ", paste0(names(args)[names(args) %in% c("returnX", "warn", "convex")], collapse = ", "), " they are set to FALSE in cv.ncvreg"))
  }
  
  ## Note ignoring ncvreg arguments
  if ("family" %in% names(args)) {
    warning("Ignoring argument 'family', only guassian family is currently supported")
  }
  
  if ("penalty" %in% names(args)) {
    warning("Ignoring argument 'penalty', only lasso penalty is currently supported")
  }
  
  if (any(c("gamma", "alpha") %in% names(args))) {
    warning(paste0("Ignoring argument(s) ", paste0(names(args)[names(args) %in% c("gamma", "alpha")], collapse = " and "), ", not used for lasso penalty"))
  }
  
  if (time) toc()
  
  original_coefs <- NULL
  ## Will select lambda, won't estimate sigma^2 without selecting lambda
  if (missing(cv_fit)) {
    if (missing(lambda) | missing(sigma2)) {
      if (missing(lambda) & missing(sigma2)) {
        if (verbose) message("Using cross validation to select lambda and estimate variance")
      } else if (missing(lambda) & !missing(sigma2)) {
        if (missing(lambda) & verbose) message("Using cross validation to select lambda")
      } else if (!missing(lambda) & missing(sigma2) & verbose) {
        message("Using cross validation to estimate variance at supplied value of lambda using linear interpolation")
      }
      
      if (time) tic(msg = "Cross Validation")
      cv.args$penalty <- "lasso"
      cv.args$X <- X
      cv.args$y <- y
      if (!(missing(lambda))) {
        lambda_max <- max(apply(ncvreg::std(X), 2, find_thresh, y))
        lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
        nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
        if (lambda_min > lambda_max | lambda > lambda_max) {
          lambda_max <- lambda + lambda / 100
          nlambda <- 2
        }
        lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
        cv.args$lambda <- lambda_seq 
      }
      if (!missing(cluster)) cv.args$cluster <- cluster ## NEED TO UPDATE
      cv_fit <- do.call("cv.ncvreg", c(cv.args, ncvreg.args))
      
      if (missing(lambda) & missing(sigma2)) {
        lambda <- cv_fit$lambda.min 
        sigma2 <- cv_fit$cve[cv_fit$min]
      } else if (missing(lambda) & !missing(sigma2)) {
        lambda <- cv_fit$lambda.min
      } else if (!missing(lambda) & missing(sigma2)) {
        if (max(cv_fit$lambda) < lambda | min(cv_fit$lambda) > lambda) stop("Supplied lambda value is outside the range of the model fit.")
        ## Make note about linear interpolation (or in documentation)
        ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
        l <- floor(ind)
        r <- ceiling(ind)
        w <- ind %% 1
        sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
      }
      
      original_coefs <- coef(cv_fit$fit, lambda = lambda)[-1]
      if (time) toc()
    } 
  } else {
    if (time) tic(msg = "Setting arguments with supplied cv.ncvreg")
    if (!missing(sigma2) & verbose) message("Overriding variance estimate in cv.ncvreg object with user specified value for sigma2.")
    if (!missing(lambda) & verbose) message("Overriding selected lambda in cv.ncvreg object with user specified value for lambda.")
    if (!missing(lambda) & missing(sigma2)) warning("Estimating variance using CV but using user specified value of lambda. Estimated variance corresponds to interpolated CVE for supplied value of lambda.")
    if (missing(lambda)) lambda <- cv_fit$lambda.min
    if (missing(sigma2)) {
      if (max(cv_fit$lambda) < lambda | min(cv_fit$lambda) > lambda) stop("Supplied lambda value is outside the range of the model fit.")
      ## Make note about linear interpolation (or in documentation)
      ind <- stats::approx(cv_fit$lambda, seq(cv_fit$lambda), lambda)$y
      l <- floor(ind)
      r <- ceiling(ind)
      w <- ind %% 1
      sigma2 <- (1-w)*cv_fit$cve[l] + w*cv_fit$cve[r]
    }
    X <- cv_fit$fit$X
    y <- cv_fit$fit$y
    original_coefs <- coef(cv_fit$fit, lambda = lambda)[-1]
    if (time) toc()
  }
  
  if (is.null(original_coefs)) {
    if (time) tic(msg = "Getting posterior mode")
    coef.args <- ncvreg.args
    coef.args$X <- X
    coef.args$y <- y
    coef.args$penalty <- "lasso"
    fit <- do.call("ncvreg", coef.args)
    original_coefs <- coef(fit, lambda = lambda)[-1]
    if (time) toc()
  }
  
  if (time) tic(msg = "Bootstrapping")  
  modes <- matrix(nrow = nboot, ncol = ncol(X))
  per_draw <- ifelse(is.character(quantiles) && quantiles == "sample", 10, length(quantiles))
  draws <- matrix(nrow = nboot * per_draw, ncol = ncol(X))
  
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "lambda", "sigma2", "significance_level", "ncvreg.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    results <- parallel::parLapply(cl=cluster, X=1:nboot, fun=bootf.r, XX=X, y=y, lambda = lambda, sigma2 = sigma2, significance_level = significance_level, ncvreg.args=ncvreg.args, rescale_original = rescale_original, quantiles = quantiles)
  }
  
  for (i in 1:nboot) {
    if (!missing(cluster)) {
      res <- results[[i]]
    } else {
      res <- bootf.r(XX=X, y=y, lambda = lambda, sigma2 = sigma2, significance_level = significance_level, ncvreg.args=ncvreg.args, rescale_original = rescale_original, quantiles = quantiles)
    }
    # print((1 + i*per_draw - per_draw):(i*per_draw))
    draws[(1 + i*per_draw - per_draw):(i*per_draw),] <- res$draws
    modes[i,] <- res$modes
  }
  if (time) toc()
  
  val <- list(draws = draws, modes = modes, estimates = original_coefs, lamdba = lambda, sigma2 = sigma2)
  
  if (returnCV) val$cv.ncvreg <- cv_fit
  
  if (time) toc() ## For overall
  structure(val, class="boot.ncvreg")
  
}
bootf.r <- function(XX, y, lambda, sigma2, significance_level = .8, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "sample") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }
  
  if (is.character(quantiles) && quantiles == "sample") {
    ps <- runif(10)  
  } else {
    ps <- quantiles
  }
  
  p <- ncol(XX)
  n <- length(y)
  
  modes <- numeric(p)
  if (time) toc()
  
  if (time) tic(msg = "Sample")
  idx_new <- sample(1:n, replace = TRUE)
  ynew <- y[idx_new]
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX[idx_new,,drop=FALSE])
  if (time) toc()
  
  if (time) tic(msg = "Lambda Sequence")
  lambda_max <- max(apply(xnew, 2, find_thresh, ynew))
  lambda_min <- lambda - lambda / 100 ## set min to be slightly smaller
  if (lambda_min > lambda_max | lambda > lambda_max) {
    lambda_max <- lambda + lambda / 100
    nlambda <- 2
  }
  ## Could use better logic to speed up
  nlambda <- ifelse(!is.null(ncvreg.args$nlambda), ncvreg.args$nlambda, 100)
  lambda_seq <- 10^(seq(log(lambda_max, 10), log(lambda_min, 10), length.out = nlambda))
  if (time) toc()
  
  if (time) tic(msg = "Fit ncvreg")
  ncvreg.args$X <- xnew
  ncvreg.args$y <- ynew
  ncvreg.args$penalty <- "lasso"
  ncvreg.args$lambda <- lambda_seq
  
  ## Ignores user specified lambda.min and nlambda
  # fit <- do.call("ncvreg", ncvreg.args[!(names(ncvreg.args) %in% c("lambda.min", "nlambda"))])
  fit <- ncvreg(xnew, ynew, penalty = "lasso", lambda = lambda_seq)
  
  coefs <- coef(fit, lambda = lambda)
  if (time) toc()
  
  if (time) tic(msg = "Compute Posterior")
  ns_index <- attr(xnew, "nonsingular")
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
  
  if (time) tic(msg = "Rescale")
  rescale <- (attr(xnew, "scale")[ns_index])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[ns_index])^(-1)
  } else {
    rescaleX <- 1
  }
  if (time) toc()
  full_rescale_factor <- rescale * rescaleX
  
  log_ps <- log(ps)
  log_one_minus_ps <- log(1 - ps)
  draws <- matrix(ncol = length(frac_lw_log), nrow = length(log_ps))
  for (i in 1:length(ps)) {
    draws[i,] <- ifelse(
      frac_lw_log >= log_ps[i],
      qnorm(log_ps[i] + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
      qnorm(log_one_minus_ps[i] + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
    ) * full_rescale_factor
  }
  
  # for (i in 1:length(frac_lw_log)) {
  #   draws[,i] <- ifelse(
  #     frac_lw_log[i] >= log_ps,
  #     qnorm(log_ps + obs_lw[i] - frac_lw_log[i], z[i] + lambda, se, log.p = TRUE),
  #     qnorm(log_one_minus_ps + obs_up[i] - frac_up_log[i], z[i] - lambda, se, lower.tail = FALSE, log.p = TRUE)
  #   )
  # }
  # 
  
  # 
  # # Function to apply qnorm based on condition
  # apply_qnorm <- function(frac_lw_log, frac_up_log, obs_lw, obs_up, z_val, full_rescale_factor, lambda_val, se_val, p) {
  #   log_p <- log(p)
  #   if (frac_lw_log >= log_p) {
  #     return(qnorm(log_p + obs_lw - frac_lw_log, mean = z_val + lambda_val, sd = se_val, log.p = TRUE))
  #   } else {
  #     return(qnorm(log(1 - p) + obs_up - frac_up_log, mean = z_val - lambda_val, sd = se_val, lower.tail = FALSE, log.p = TRUE))
  #   }
  # }
  # 
  # # Apply the function for each combination of ps and frac_lw_log/obs_lw/etc.
  # draws <- outer(1:length(ps), 1:length(frac_lw_log), Vectorize(function(i, j) {
  #   apply_qnorm(frac_lw_log[j], frac_up_log[j], obs_lw[j], obs_up[j], z[j], full_rescale_factor[j], lambda, se, ps[i])
  # }))

  if (time) toc()
  
  if (time) tic(msg = "Return result")
  
  # draws <- sapply(1:ncol(draws), function(x) draws[,x]*full_rescale_factor[x])
  # draws <- draws*matrix(full_rescale_factor, ncol = length(frac_lw_log), nrow = length(ps), byrow = TRUE)
  modes[ns_index] <- (modes * rescale) * rescaleX
  if (length(ns_index) < ncol(draws)) draws[,!(1:ncol(draws) %in% ns_index)] <- rep(NA, length(ps))
  modes[!(1:length(modes) %in% ns_index)] <- NA
  
  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  if (time) toc()
  if (time) toc()
  return(ret)
  
}
find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }


