# 'number': if you want a summary for more than S variables
local_mfdr <- function(fit, lambda, number=NULL, cutoff=NULL, X=NULL, y=NULL) {

  # Check for valid args
  if(!is.null(number) && (number < 1)) stop("'number' should be a positive integer")
  if(!is.null(cutoff) && (cutoff > 1 | cutoff <= 0)) stop("'cutoff' should be in the interval (0,1]")

  # Extract standardized X, y
  if (is.null(X) & is.null(fit$X)) {
      stop("This procedure requires X and y.  Either supply X and y, or fit the model using the option 'returnX = TRUE'")
  }
  if(class(fit)[1] == "ncvsurv") {
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

  if(class(fit)[1] == "ncvreg") {
    if(fit$family == "gaussian") {
      # Linear Regression
      bb <- beta[-1][ns]*sc
      r <- yy - XX %*% bb
      z <- crossprod(XX, r)/n + bb
      rss <- approxfun(fit$lambda, fit$loss)
      sig.est <- sqrt(rss(lambda)/(n - S + 1))
      z <- z/(sig.est/sqrt(n))
    } else if (fit$family == "binomial") {
      # Logistic regression
      bb <- beta[-1][ns]*sc
      a <- sum(cn*beta[-1][ns]) + beta[1]

      # Setup the score vector and W matrix
      P <- 1/(1 + exp(-a-(XX %*% bb)))
      U <- yy - P
      W <- diag(as.vector(P*(1 - P)))

      # Calculate v_j and z_j
      z <- numeric(p)
      for (j in 1:p){
        vj <- t(XX[,j]) %*% W %*% XX[,j]
        z[j] <- (XX[,j] %*% U + vj * bb[j])/(sqrt(vj))  ### j+1 bc intercept
      }
    }
  }

  # Cox regression
  if(class(fit)[1] == "ncvsurv"){
    # Setup standardized beta
    bb <- beta[ns]*sc

    # Calculate score vector and W (maybe diagonalize it for speed?)
    ind <- approx(fit$lambda, seq(fit$lambda), lambda)$y
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
    z <- numeric(p)
    for (j in 1:p){
      vj <- t(XX[,j]) %*% W %*% XX[,j]
      z[j] <- (XX[,j] %*% U + vj * bb[j])/(sqrt(vj))
    }
  }

  # Calculate locfdr
  f <- density(z[pen.idx])
  ff <- approxfun(f$x, f$y)
  est.gam <- pmin(dnorm(z[pen.idx], 0, 1)/ff(z[pen.idx]), 1)

  # setup results
  Estimate <- if (class(fit)[1] == "ncvreg") beta[-1][ns][pen.idx] else beta[ns][pen.idx]
  results <- data.frame(Estimate = Estimate, z = z[pen.idx], mfdr = est.gam)
  rownames(results) <- names(Estimate)

  # Results for unpenalized vars
  # This uses the effect of the high dim selections as an offset
  unpen.res <- NULL
  if(sum(pen.idx) < length(pen.idx)){
    if(class(fit)[1] == "ncvreg"){
      off <- XX[,pen.idx] %*% bb[pen.idx]
      unpen.res <- summary(glm(yy ~ XX[,!pen.idx], offset = off, family = fit$family))$coef
      unpen.res <- data.frame(Estimate = beta[-1][ns][!pen.idx], std.error = unpen.res[-1,2]/sc[ns][!pen.idx], statistic = unpen.res[-1,3], p.value = unpen.res[-1,4])
      rownames(unpen.res) <- names(bb)[!pen.idx]
    } else {
      off <- XX[,pen.idx] %*% bb[pen.idx]
      unpen.res <- summary(survival::coxph(survival::Surv(y,d) ~ XX[,!pen.idx] + offset(off)))$coefficients
      unpen.res <- data.frame(Estimate = unpen.res[,1]/sc[ns][!pen.idx], std.error = unpen.res[,3]/sc[ns][!pen.idx], statistic = unpen.res[,4], p.value = unpen.res[,5])
      rownames(unpen.res) <- names(bb)[!pen.idx]
    }
  }

  # If number and cutoff are both unspecified return results for selected vars
  results <- results[order(results$mfdr),]
  if(is.null(number) & is.null(cutoff)){
    results <- results[results$Estimate != 0,]
    return(list(pen.vars = results[!is.na(results$mfdr),],unpen.vars = unpen.res ))
  }

  # If cutoff is null return 1:number, else return the minimum
  if(is.null(cutoff)) {
    nShow <- min(number, p)
  } else {
    nShow <- min(sum(results$mfdr <= cutoff, na.rm = TRUE), number, p)
  }
  results$Selected <- ifelse(results$Estimate != 0, "*"," ")
  return(list(pen.vars=results[1:nShow,], unpen.vars=unpen.res))
}
