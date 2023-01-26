#' Estimate local mFDR for all features
#' 
#' `local_mfdr()` is called by `summary.ncvreg()`, which typically offers a more convenient interface to users.
#' If, however, you are working with local mfdrs programmatically rather than interactively, you probably want to
#' use `local_mfdr()`, which skips the sorting, filtering, and print formatting of `summary.ncvreg()`.
#' 
#' @param fit      A fitted `ncvreg` or `ncvsurv` object.
#' @param lambda   The value of lambda at which inference should be carried out.
#' @param X,y      The design matrix and response used to fit the model; in most cases, it is not necessary to provide
#'   `X` and `y` as they are returned by `ncvreg`, but see the `returnX` argument in [ncvreg()].
#' @param method   What method should be used to calculate the local fdr?  Options are `ashr` (which tends to be more
#'   accurate) and `kernel` (which requires no additional packages).  The default is to use `ashr` if the package is
#'   installed.
#' @param sigma    For linear regression models, users can supply an estimate of the residual standard deviation.
#'   The default is to use RSS / DF, where degrees of freedom are approximated using the number of nonzero coefficients.
#' @param ...      Additional arguments to `ash()` if using `method='ashr'`.
#'   
#' @return If all features are penalized, then the object returns a data frame with one row per feature and four columns:
#' * `Estimate`: The coefficient estimate from the penalized regression fit
#' * `z`: A test statistic that approximately follows a standard normal distribution under the null hypothesis that the
#'        feature is marginally independent of the outcome
#' * `mfdr`: The estimated marginal local false discovery rate
#' * `Selected`: Features with nonzero coefficient estimates are given an asterisk
#' 
#' If some features are penalized and others are not, then a list is returned with two elements: `pen.vars`, which consists
#' of the data frame described above, and `unpen.vars`, a data frame with four columns: `Estimate`, `SE`, `Statistic`, and
#' `p.value`.  The standard errors and p-values are based on a classical `lm`/`glm`/`coxph` model using the effect of the
#' penalized features as an offset.
#' 
#' @seealso [summary.ncvreg()]
#'
#' @examples
#' # Linear regression
#' data(Prostate)
#' fit <- ncvreg(Prostate$X, Prostate$y)
#' local_mfdr(fit, 0.1)
#' 
#' fit <- ncvreg(Prostate$X, Prostate$y, penalty.factor=rep(0:1, each=4))
#' local_mfdr(fit, 0.1)
#' 
#' # Logistic regression
#' data(Heart)
#' X <- Heart$X
#' y <- Heart$y
#' fit <- ncvreg(X, y, family='binomial')
#' local_mfdr(fit, 0.1)
#' 
#' # Cox regression
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#' fit <- ncvsurv(X, y)
#' local_mfdr(fit, 0.1)
#' @export

local_mfdr <- function(fit, lambda, X=NULL, y=NULL, method=c('ashr', 'kernel'), sigma, ...) {
  
  # Determine method, if missing
  if (!inherits(fit, 'ncvreg')) stop('"fit" must be an ncvreg or ncvsurv object', call.=FALSE)
  if (missing(method)) {
    if (requireNamespace('ashr', quietly=TRUE)) {
      method <- 'ashr'
    } else {
      method <- 'kernel'
      if (is.null(getOption('ncvreg.ashr.warn'))) {
        message('Using a basic kernel estimate for local fdr; consider installing the ashr package for more accurate estimation.  See ?local_mfdr')
        options(ncvreg.ashr.warn = FALSE)
      }
    }
  }

  # Extract standardized X, y
  if (is.null(X) & is.null(fit$X)) {
      stop("This procedure requires X and y.  Either supply X and y, or fit the model using the option 'returnX = TRUE'", call.=FALSE)
  }
  if (inherits(fit, "ncvsurv")) {
    tmp <- if (is.null(fit$X)) ncvsurv(X, y) else fit
    XX <- tmp$X
    y <- tmp$time
    d <- tmp$fail
  } else {
    tmp <- if (is.null(fit$X)) ncvreg(X, y, family=fit$family) else fit
    XX <- tmp$X
    yy <- tmp$y
  }

  # Setup general
  ns <- attr(XX, "nonsingular")
  sc <- attr(XX, "scale")[ns]
  cn <- attr(XX, "center")[ns]
  n <- nrow(XX)
  p <- ncol(XX)
  S <- predict(fit, type = "nvars", lambda = lambda)
  beta <- coef(fit, lambda=lambda)
  pen.idx <- fit$penalty.factor > 0

  if (!inherits(fit, "ncvsurv")) {
    if (fit$family == "gaussian") {
      # Linear Regression
      bb <- beta[-1][ns]*sc
      r <- yy - XX %*% bb
      z <- crossprod(XX, r)/n + bb
      if (missing(sigma)) {
        rss <- stats::approxfun(fit$lambda, fit$loss)
        sigma <- sqrt(rss(lambda)/(n - S + 1))
        if (S > n) stop('Default estimate of sigma, sqrt(RSS/DF), is invalid because DF is negative. Supply an estimate or use cross-validation.', call.=FALSE)
      }
      z <- z/(sigma/sqrt(n))
    } else if (fit$family == "binomial") {
      # Logistic regression
      bb <- beta[-1][ns]*sc
      a <- sum(cn*beta[-1][ns]) + beta[1]

      # Setup the score vector and W matrix
      P <- 1/(1 + exp(-a-(XX %*% bb)))
      U <- yy - P
      W <- diag(as.vector(P*(1 - P)))

      # Calculate v_j and z_j
      z <- double(p)
      for (j in 1:p){
        vj <- t(XX[,j]) %*% W %*% XX[,j]
        z[j] <- (XX[,j] %*% U + vj * bb[j])/(sqrt(vj))  ### j+1 bc intercept
      }
    } else if (fit$family == "poisson") {
      # Poisson regression
      bb <- beta[-1][ns]*sc
      a <- sum(cn*beta[-1][ns]) + beta[1]
      
      # Setup the score vector and W matrix
      mu <- exp(a+(XX %*% bb))
      U <- yy - mu
      W <- diag(as.vector(mu))
      
      # Calculate v_j and z_j
      z <- double(p)
      for (j in 1:p){
        vj <- t(XX[,j]) %*% W %*% XX[,j]
        z[j] <- (XX[,j] %*% U + vj * bb[j])/(sqrt(vj))  ### j+1 bc intercept
      }
    }
  }
  
  # Cox regression
  if (inherits(fit, "ncvsurv")) {
    # Setup standardized beta
    bb <- beta[ns]*sc

    # Calculate score vector and W (maybe diagonalize it for speed?)
    ind <- stats::approx(fit$lambda, seq(fit$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    x <- ind %% 1
    Eta <- (1-x)*fit$Eta[,l] + x*fit$Eta[,r]
    rsk <- rev(cumsum(rev(exp(Eta))))
    P <- outer(exp(Eta), rsk, '/')
    P[upper.tri(P)] <- 0
    U <- d - P%*%d
    W <- -P %*% diag(d) %*% t(P)
    diag(W) <- diag(P %*% diag(d) %*% t(1-P))

    # Calculate v_j and z_j
    z <- double(p)
    for (j in 1:p){
      vj <- t(XX[,j]) %*% W %*% XX[,j]
      z[j] <- (XX[,j] %*% U + vj * bb[j])/(sqrt(vj))
    }
  }

  # Calculate locfdr
  if (method=='ashr') {
    ash_fit <- ashr::ash(z[pen.idx], rep(1, sum(pen.idx)), optmethod='mixEM', ...)
    est.gam <- ashr::get_lfdr(ash_fit)
  } else {
    f <- density(z[pen.idx])
    ff <- stats::approxfun(f$x, f$y)
    est.gam <- pmin(dnorm(z[pen.idx], 0, 1)/ff(z[pen.idx]), 1)
  }    

  # Setup results and return, if no unpenalized variables
  Estimate <- if (inherits(fit, "ncvsurv")) beta[ns][pen.idx] else beta[-1][ns][pen.idx]
  results <- data.frame(Estimate = Estimate, z = z[pen.idx], mfdr = est.gam)
  rownames(results) <- names(Estimate)
  results$Selected <- ifelse(results$Estimate != 0, "*"," ")
  if (sum(pen.idx) == length(pen.idx)) return(results)

  # Results for unpenalized vars, using the effect of the high dim selections as an offset
  if (!inherits(fit, "ncvsurv")) {
    off <- XX[, pen.idx] %*% bb[pen.idx]
    unpen.res <- summary(glm(yy ~ XX[,!pen.idx], offset = off, family = fit$family))$coef
    unpen.res <- data.frame(Estimate = beta[-1][ns][!pen.idx], std.error = unpen.res[-1,2]/sc[ns][!pen.idx], statistic = unpen.res[-1,3], p.value = unpen.res[-1,4])
    rownames(unpen.res) <- names(bb)[!pen.idx]
  } else {
    off <- XX[, pen.idx] %*% bb[pen.idx]
    unpen.res <- summary(survival::coxph(survival::Surv(y,d) ~ XX[,!pen.idx] + offset(off)))$coefficients
    unpen.res <- data.frame(Estimate = unpen.res[,1]/sc[ns][!pen.idx], std.error = unpen.res[,3]/sc[ns][!pen.idx], statistic = unpen.res[,4], p.value = unpen.res[,5])
    rownames(unpen.res) <- names(bb)[!pen.idx]
  }
  list(pen.vars=results, unpen.vars=unpen.res)
}
