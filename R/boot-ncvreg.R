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
#' @export boot.ncvreg
boot.ncvreg <- function(X, y, cv_fit, lambda, sigma2, nboot = 100, ..., cluster, seed, returnCV=FALSE, verbose = TRUE, time = FALSE, quantiles = "sample") {
  
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
  per_draw <- ifelse(is.character(quantiles), ifelse(quantiles == "sample", 1, 1), length(quantiles))
  draws <- matrix(nrow = nboot * per_draw, ncol = ncol(X))
  
  if (is.character(quantiles) && quantiles == "sample" | !is.character(quantiles)) {
    bootf <- bootf1 
  } else if (is.character(quantiles) && quantiles == "disturbed") {
    bootf <- bootf2
  } else if (is.character(quantiles) && quantiles == "mode") {
    bootf <- bootf3
  } else if (is.character(quantiles) && quantiles == "zs") {
    bootf <- bootf4
  } else if (is.character(quantiles) && quantiles == "zerosample") {
    bootf <- bootf5
  }
  
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "lambda", "sigma2", "ncvreg.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    results <- parallel::parLapply(cl=cluster, X=1:nboot, fun=bootf, XX=X, y=y, lambda = lambda, sigma2 = sigma2, ncvreg.args=ncvreg.args, rescale_original = rescale_original, quantiles = quantiles)
  }
  
  for (i in 1:nboot) {
    if (!missing(cluster)) {
      res <- results[[i]]
    } else {
      res <- bootf(XX=X, y=y, lambda = lambda, sigma2 = sigma2, ncvreg.args=ncvreg.args, rescale_original = rescale_original, quantiles = quantiles)
    }
    draws[(1 + i*per_draw - per_draw):(i*per_draw),] <- res$draws
    modes[i,] <- res$modes
  }
  if (time) toc()
  
  val <- list(draws = draws, modes = modes, estimates = original_coefs, lambda = lambda, sigma2 = sigma2)
  
  if (returnCV) val$cv.ncvreg <- cv_fit
  
  if (time) toc() ## For overall
  structure(val, class="boot.ncvreg")
  
}
bootf1 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "sample") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }
  
  if (is.character(quantiles) && quantiles == "sample") {
    ps <- runif(1)  
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
  nonsingular <- attr(xnew, "nonsingular")
  
  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  if (time) toc()
  full_rescale_factor <- rescale * rescaleX
  
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
  
  log_ps <- log(ps) 
  log_one_minus_ps <- log(1 - ps)
  draws <- matrix(ncol = p, nrow = length(log_ps))
  for (i in 1:length(ps)) {
    draws[i,nonsingular] <- ifelse(
      frac_lw_log >= log_ps[i],
      qnorm(log_ps[i] + obs_lw - frac_lw_log, z + lambda, se, log.p = TRUE),
      qnorm(log_one_minus_ps[i] + obs_up - frac_up_log, z - lambda, se, lower.tail = FALSE, log.p = TRUE)
    ) * full_rescale_factor
  }

  if (time) toc()
  
  if (time) tic(msg = "Return result")
  
  modes[nonsingular] <- (modes * rescale) * rescaleX
  if (length(nonsingular) < ncol(draws)) draws[,!(1:ncol(draws) %in% nonsingular)] <- rep(NA, length(ps))
  modes[!(1:length(modes) %in% nonsingular)] <- NA
  
  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  if (time) toc()
  if (time) toc()
  return(ret)
  
}
bootf2 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "disturbed") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
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
  nonsingular <- attr(xnew, "nonsingular")
  
  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  if (time) toc()
  full_rescale_factor <- rescale * rescaleX
  
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
  
  if (time) toc()
  
  if (time) tic(msg = "Return result")
  
  modes[nonsingular] <- (modes * rescale) * rescaleX
  
  if (length(nonsingular) < ncol(draws)) draws[,!(1:ncol(draws) %in% nonsingular)] <- NA
  modes[!(1:length(modes) %in% nonsingular)] <- NA
  
  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  if (time) toc()
  if (time) toc()
  return(ret)
  
}
bootf3 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "mode") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
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
  nonsingular <- attr(xnew, "nonsingular")
  
  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  if (time) toc()
  full_rescale_factor <- rescale * rescaleX
  
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

  modes <- coefs[-1] ## Coefs only returned for nonsingular columns of X
  modes[nonsingular] <- (modes * rescale) * rescaleX
  modes[!(1:length(modes) %in% nonsingular)] <- NA
  
  ret <- list(modes, modes)
  names(ret) <- c("draws", "modes")

  return(ret)
  
}
bootf4 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "zs") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
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
  nonsingular <- attr(xnew, "nonsingular")
  
  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1)
  } else {
    rescaleX <- 1
  }
  if (time) toc()
  full_rescale_factor <- rescale * rescaleX
  
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
bootf5 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "zerosample") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }

  p <- ncol(XX)
  n <- length(y)
  
  if (time) toc()
  
  if (time) tic(msg = "Sample")
  idx_new <- sample(1:n, replace = TRUE)
  ynew <- y[idx_new]
  ynew <- ynew - mean(ynew)
  xnew <- ncvreg::std(XX[idx_new,,drop=FALSE])
  nonsingular <- attr(xnew, "nonsingular")
  
  rescale <- (attr(xnew, "scale")[nonsingular])^(-1)
  if (!is.null(attr(XX, "scale")) & rescale_original) {
    rescaleX <-  (attr(XX, "scale")[nonsingular])^(-1) ## may need to take a closer look at this
  } else {
    rescaleX <- 1
  }
  if (time) toc()
  full_rescale_factor <- rescale * rescaleX
  
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

  if (time) toc()
  
  if (time) tic(msg = "Return result")
  
  tmp <- modes
  modes <- numeric(p)
  modes[nonsingular] <- tmp * full_rescale_factor 
  if (length(nonsingular) < ncol(draws)) draws[1,!(1:ncol(draws) %in% nonsingular)] <- NA
  modes[!(1:length(modes) %in% nonsingular)] <- NA
  
  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  if (time) toc()
  if (time) toc()
  return(ret)
  
}
bootf_nosample <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "disturbed") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }
  
  p <- ncol(XX)
  n <- length(y)
  
  modes <- numeric(p)
  if (time) toc()
  
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
  
  if (time) toc()
  
  if (time) tic(msg = "Return result")
  
  modes[nonsingular] <- (modes * rescale) * rescaleX
  if (length(nonsingular) < ncol(draws)) draws[,!(1:ncol(draws) %in% nonsingular)] <- NA
  modes[!(1:length(modes) %in% nonsingular)] <- NA
  
  ret <- list(draws, modes)
  names(ret) <- c("draws", "modes")
  if (time) toc()
  if (time) toc()
  return(ret)
  
}
bootf_nosample2 <- function(XX, y, lambda, sigma2, ncvreg.args, rescale_original = TRUE, time = FALSE, quantiles = "disturbed") {
  
  if (time) tic(msg = "Overall - Bootstrap")
  if (time) tic(msg = "Prep")
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }
  
  p <- ncol(XX)
  n <- length(y)
  
  modes <- numeric(p)
  if (time) toc()
  
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
find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }


