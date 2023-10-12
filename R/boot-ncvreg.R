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
#' @param cvncvreg Instead of \code{X} and \code{y}, an object of type 
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
#' @param trace If set to TRUE, inform the user of progress by displaying 
#' progress through bootstrapping and by announcing the beginning of each CV 
#' fold if \code{cv.ncvreg} is ran.  Default is FALSE.
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
#' ## Specifying cvncvreg
#' tmp <- boot.ncvreg(cvncvreg = cv.ncvreg(dat$X, dat$y, penalty = "lasso", returnX = TRUE))
#' 
#' @export boot.ncvreg
boot.ncvreg <- function(X, y, cvncvreg, lambda, sigma2, significance_level = 0.8, nboot = 100, ..., cluster, seed, returnCV=FALSE, trace=FALSE, verbose = TRUE) {
  
  ## name conflict: lambda -> removed from being able to specify own sequence fo cv.ncvreg as temp solution
  ## Handle setting 2 seeds or resetting seeds? Don't really want to support this -> set seed one
  ## Both can be circumvented by using cv.ncvreg route if desired
  ## Best naming convention for cvncvreg?
  ## Should warn / convex be saved by default?
  ## Match arguments (so they don't need to be named??)
  ## NEED TO CHECK IF MISSING ONLY BEFORE I ACTIVELY ASSIGN SOMETHING
  
  if ((missing(X) | missing(y)) & (missing(cvncvreg) || class(cvncvreg) != "cv.ncvreg")) {
    stop("Either X and y or an object of class cv.ncvreg must be supplied.")
  }
  
  if (!missing(cvncvreg)) {
    if (!all(c("X", "y") %in% names(cvncvreg$fit))) {
      stop("fit object in cvncvreg missing X and y, please rerun cv.ncvreg with returnX = TRUE or supply X and y directly (without also specifying a cv.ncvreg object)") 
    }
    
    if(cvncvreg$fit$penalty != "lasso") {
      stop(paste0("cv.ncvreg fit with ", cvncvreg$fit$penalty, " penalty, but only 'lasso' penalty is currently supported for boot.ncvreg"))
    }
    
    if(cvncvreg$fit$family != "gaussian") {
      stop(paste0("cv.ncvreg fit with ", cvncvreg$fit$family, " family, but only 'gaussian' family is currently supported for boot.ncvreg"))
    }
  }
  
  if (!missing(y) & !missing(cvncvreg)) {
    warning("Ignoring supplied values of y, using y from supplied cv.ncvreg object")
  }
  
  if (!missing(X) & !missing(cvncvreg)) {
    if ("X" %in% names(cvncvreg$fit)) {
      warning("Ignoring supplied values of X, using X from supplied cv.ncvreg object") 
    }
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
  if (length(args) > 0 & !missing(cvncvreg)) {
    warning("Additional arguments are ignored when cv.ncvreg object supplied")
  }
  
  if (length(args) > 0 & any(is.null(names(args)))) {
    stop("Please supply names for all additional arguments passed to ...")
  }
  
  
  cv.args <- args[names(args) %in% c("nfolds", "fold", "returnY")]
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
  
  original_coefs <- NULL
  ## Will select lambda, won't estimate sigma^2 without selecting lambda
  if (missing(cvncvreg)) {
    if (missing(lambda) | missing(sigma2)) {
      if (missing(lambda) & missing(sigma2)) {
        if (verbose) message("Using cross validation to select lambda and estimate variance")
      } else if (missing(lambda) & !missing(sigma2)) {
        if (missing(lambda) & verbose) message("Using cross validation to select lambda")
      } else if (!missing(lambda) & missing(sigma2) & verbose) {
        message("Using cross validation to estimate variance at supplied value of lambda using linear interpolation")
      }
      
      cv.args$trace <- trace
      cv.args$penalty <- "lasso"
      cv.args$X <- X
      cv.args$y <- y
      if (!missing(cluster)) cv.args <- cluster
      cvncvreg <- do.call("cv.ncvreg", c(cv.args, ncvreg.args))
      if (missing(lambda)) lambda <- cvncvreg$lambda.min
      if (missing(sigma2)) {
        if (max(cvncvreg$lambda) < lambda | min(cvncvreg$lambda) > lambda) stop("Supplied lambda value is outside the range of the model fit.")
        ## Make note about linear interpolation (or in documentation)
        ind <- stats::approx(cvncvreg$lambda, seq(cvncvreg$lambda), lambda)$y
        l <- floor(ind)
        r <- ceiling(ind)
        w <- ind %% 1
        sigma2 <- (1-w)*cvncvreg$cve[l] + w*cvncvreg$cve[r]
      }
      original_coefs <- coef(cvncvreg$fit, lambda = lambda)[-1]
    } 
  } else {
    if (!missing(sigma2) & verbose) message("Overriding variance estimate in cv.ncvreg object with user specified value for sigma2.")
    if (!missing(lambda) & verbose) message("Overriding selected lambda in cv.ncvreg object with user specified value for lambda.")
    if (!missing(lambda) & missing(sigma2)) warning("Estimating variance using CV but using user specified value of lambda. Estimated variance corresponds to interpolated CVE for supplied value of lambda.")
    if (missing(lambda)) lambda <- cvncvreg$lambda.min
    if (missing(sigma2)) {
      if (max(cvncvreg$lambda) < lambda | min(cvncvreg$lambda) > lambda) stop("Supplied lambda value is outside the range of the model fit.")
      ## Make note about linear interpolation (or in documentation)
      ind <- stats::approx(cvncvreg$lambda, seq(cvncvreg$lambda), lambda)$y
      l <- floor(ind)
      r <- ceiling(ind)
      w <- ind %% 1
      sigma2 <- (1-w)*cvncvreg$cve[l] + w*cvncvreg$cve[r]
    }
    X <- cvncvreg$fit$X
    y <- cvncvreg$fit$y
    original_coefs <- coef(cvncvreg$fit, lambda = lambda)[-1]
  }
  
  if (is.null(original_coefs)) {
    coef.args <- ncvreg.args
    coef.args$X <- X
    coef.args$y <- y
    coef.args$penalty <- "lasso"
    fit <- do.call("ncvreg", coef.args)
    original_coefs <- coef(fit, lambda = lambda)[-1]
  }
  
  lowers <- uppers <- modes <- matrix(nrow = nboot, ncol = ncol(X))

  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "lambda", "sigma2", "significance_level", "ncvreg.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    results <- parallel::parLapply(cl=cluster, X=1:nboot, fun=bootf, XX=X, y=y, lambda = lambda, sigma2 = sigma2, significance_level = significance_level, ncvreg.args=ncvreg.args)
  }
  
  for (i in 1:nboot) {
    if (!missing(cluster)) {
      res <- results[[i]]
    } else {
      res <- bootf(XX=X, y=y, lambda = lambda, sigma2 = sigma2, significance_level = significance_level, ncvreg.args=ncvreg.args)
    }
    lowers[i,] <- res$lowers
    uppers[i,] <- res$uppers
    modes[i,] <- res$modes
  }
  
  val <- list(lowers = lowers, uppers = uppers, modes = modes, estimates = original_coefs)
  
  if (returnCV) val$cv.ncvreg <- cvncvreg
  structure(val, class="boot.ncvreg")
}
bootf <- function(XX, y, lambda, sigma2, significance_level = .8, ncvreg.args) {
  
  if (missing(ncvreg.args)) {
    ncvreg.args <- list()
  }
  
  lower_p <- (1 - significance_level) / 2
  upper_p <- significance_level + lower_p
  p <- ncol(XX)
  n <- length(y)
  
  modes <- uppers <- lowers <- numeric(p)
  
  idx_new <- sample(1:n, replace = TRUE)
  ynew <- y[idx_new]
  xnew <- ncvreg::std(XX[idx_new,,drop=FALSE])
  
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
  fit <- do.call("ncvreg", ncvreg.args[!(ncvreg.args %in% c("lambda.min", "nlambda"))])
  
  coefs <- coef(fit, lambda = lambda)
  
  ns_index <- attr(xnew, "nonsingular")
  modes <- coefs[-1] ## Coefs only returned for nonsingular columns of X
  
  partial_residuals <-  ynew - (coefs[1] + as.numeric(xnew %*% modes) - (xnew * matrix(modes, nrow=nrow(xnew), ncol=ncol(xnew), byrow=TRUE)))
  
  z <- (1/n)*colSums(xnew * partial_residuals)
  se <- sqrt(sigma2 / n)
  
  obs_lw <- pnorm(0, z + lambda, se)
  obs_up <- pnorm(0, z - lambda, se, lower.tail = FALSE)
  
  # C <- exp(n*z*lambda / sigma2)^2
  # lower_adj <- obs_lw + obs_up*(C^(-1))
  # upper_adj <- obs_up + obs_lw*C
  
  C <- exp(n*z*lambda / sigma2)
  t <- obs_lw * C + obs_up * C^(-1)
  lower_adj <- t / C 
  upper_adj <- t / (C^-1)
  
  prop_lw <- obs_lw  / lower_adj
  
  lower <- ifelse(
    prop_lw >= lower_p,
    qnorm(lower_p * lower_adj, z + lambda, se),
    qnorm(upper_p * upper_adj, z - lambda, se, lower.tail = FALSE)
  )
  upper <- ifelse(
    prop_lw >= upper_p,
    qnorm(upper_p * lower_adj, z + lambda, se),
    qnorm(lower_p * upper_adj, z - lambda, se, lower.tail = FALSE)
  )
  
  rescale <- (attr(xnew, "scale")[ns_index])^(-1)
  if (!is.null(attr(XX, "scale"))) {
    rescaleX <-  (attr(XX, "scale")[ns_index])^(-1)
  } else {
    rescaleX <- 1
  }
  
  lowers[ns_index] <- (lower * rescale) * rescaleX
  uppers[ns_index] <- (upper * rescale) * rescaleX
  modes[ns_index] <- (modes * rescale) * rescaleX
  lowers[!(1:length(lowers) %in% ns_index)] <- NA
  uppers[!(1:length(uppers) %in% ns_index)] <- NA
  modes[!(1:length(modes) %in% ns_index)] <- NA
  
  ret <- list(lowers, uppers, modes)
  names(ret) <- c("lowers", "uppers", "modes")
  return(ret)
  
}
find_thresh <- function(x, y) { abs(t(x) %*% y) / length(y) }


