# 'number': if you want a summary for more than S variables
local_mfdr <- function(fit, lambda, X = NULL, y = NULL, number=NULL, cutoff=NULL){
  
  # Check for valid args
  if(!is.null(number)){
    if(number < 1) stop("'number' should be a positive integer")
  }
  if(!is.null(cutoff)){
    if(cutoff > 1 | cutoff <= 0) stop("'cutoff' should be in the interval (0,1]")
  }  
  if(class(fit)[1] == "ncvreg_raw") {
    stop("local_mfdr currently does not support raw data fit with ncvreg_raw")
  }
  # Standardize X
  if(is.null(X)) {
    if(is.null(fit$X)) {
      stop("This procedure requires X and y, either supply X and y, or fit the model using the option 'returnX = TRUE'")
    } else {
      XX <- fit$X
      y <- fit$y
    }
  } else {
    if(class(fit)[1] == "ncvsurv") { 
      XX <- std(X[fit$order,])
    } else {
      XX <- std(X)
    }
  }
  
  ### Setup general
  ns <- attr(XX, "nonsingular")
  sc <- attr(XX, "scale")
  cen <- attr(XX, "center")
  n <- nrow(XX)
  p <- ncol(XX)
  S <- predict(fit, type = "nvars", lambda = lambda)
  beta <- coef(fit, lambda=lambda)
  pen.idx <- fit$penalty.factor > 0
  
  ##### Linear Regression
  if(class(fit)[1] == "ncvreg"){
    if(fit$family == "gaussian"){
      ### Setup standardized beta and centered y
      yy <- y - mean(y)
      bb <- c(mean(y), beta[ns+1]*sc)
      
      ### Calculate standardized z_j's
      R <- yy - cbind(1, XX) %*% bb
      z <- (1/n)*t(XX) %*% R + bb[-1]
      rss <- approxfun(fit$lambda, fit$loss)
      sig.est <- sqrt(rss(lambda)/(n - S + 1))
      z <- z/(sig.est/sqrt(n))
    }
    
    ### Logistic regression 
    else if (fit$family == "binomial"){
      ## Setup standardized beta
      yy <- y
      if(is.factor(y)){
        yy <- as.numeric(y)-1
      } 
      bb <- c(beta[1] + sum(beta[ns+1]*cen), beta[ns+1]*sc)
      
      ### Setup the score vector and W matrix
      P <- 1/(1 + exp(-cbind(1,XX) %*% bb))
      U <- yy - P   
      W <- diag(as.vector(P*(1 - P)))  
      
      ### Calculate v_j and z_j
      z <- numeric(p)
      for (j in 1:p){
        vj <- t(XX[,j]) %*% W %*% XX[,j]
        z[j] <- (XX[,j] %*% U + vj * bb[j+1])/(sqrt(vj))  ### j+1 bc intercept
      }
    }
  }
  
  ### Cox regression
  if(class(fit)[1] == "ncvsurv"){
    ## Setup standardized beta
    bb <- beta[ns]*sc
    
    ### Calculate score vector and W (maybe diagonalize it for speed?)
    d <- fit$fail
    rsk <- rev(cumsum(rev(exp(fit$Eta[,lid]))))
    P <- outer(exp(fit$Eta[,lid]), rsk, '/')
    P[upper.tri(P)] <- 0
    U <- d - P%*%d  
    W <- -P %*% diag(d) %*% t(P)
    #W <- matrix(0,n,n)
    diag(W) <- diag(P %*% diag(d) %*% t(1-P))
    
    ### Calculate v_j and z_j
    z <- numeric(p)
    for (j in 1:p){
      vj <- t(XX[,j]) %*% W %*% XX[,j]
      z[j] <- (XX[,j] %*% U + vj * bb[j])/(sqrt(vj))
    }
  }
  
  ### Calculate locfdr
  if (require(ashr)) {
    ash_fit <- ash(z[pen.idx], rep(1, sum(pen.idx)), optmethod='mixEM')
    est.gam <- get_lfdr(ash_fit)
  } else {
    f <- density(z[pen.idx])
    ff <- approxfun(f$x, f$y)
    est.gam <- pmin(dnorm(z[pen.idx], 0, 1)/ff(z[pen.idx]), 1)
  }
  
  #### Calculate Fdr (using both est cdf and empirical cdf)
  #### Remove this section if only care about locfdr
  # est.Fdr <- emp.Fdr <- numeric(p)
  # est.Fdr[!pen.idx] <- emp.Fdr[!pen.idx] <- NA
  # for(j in (1:p)[pen.idx]){
  #   if(z[j] < 0) {
  #     emp.Fdr[j] <- pnorm(z[j])/(sum(z <= z[j])/p)
  #   } else {
  #     emp.Fdr[j] <- (1-pnorm(z[j]))/(sum(z >= z[j])/p)
  #   }
  # }
  
  ### setup results
  Estimate <- if (class(fit)[1] == "ncvreg") beta[-1][pen.idx] else beta[pen.idx]
  results <- data.frame(Estimate = Estimate, z = z[pen.idx], mfdr = est.gam)
  rownames(results) <- names(Estimate)

  ### Results for unpenalized vars
  ### This uses the effect of the high dim selections as an offset
  unpen.res <- NULL
  if(sum(pen.idx) < length(pen.idx)){
    if(class(fit)[1] == "ncvreg"){
      off <- XX[,pen.idx] %*% bb[c(FALSE, pen.idx)]
      unpen.res <- summary(glm(yy ~ XX[,!pen.idx], offset = off, family = fit$family))$coef
      unpen.res <- data.frame(Estimate = unpen.res[-1,1]/sc[!pen.idx], std.error = unpen.res[-1,2]/sc[!pen.idx], statistic = unpen.res[-1,3], p.value = unpen.res[-1,4])
      rownames(unpen.res) <- names(beta)[c(FALSE, !pen.idx)]
    } else {
      off <- XX[,pen.idx] %*% bb[pen.idx]
      unpen.res <- summary(survival::coxph(survival::Surv(y,d) ~ XX[,!pen.idx] + offset(off)))$coefficients
      unpen.res <- data.frame(Estimate = unpen.res[,1]/sc[!pen.idx], std.error = unpen.res[,3]/sc[!pen.idx], statistic = unpen.res[,4], p.value = unpen.res[,5])
      rownames(unpen.res) <- names(beta)[!pen.idx]
    }
  }
  
  ### If number and cutoff are both unspecified return results for selected vars
  results <- results[order(results$mfdr),]
  if(is.null(number) & is.null(cutoff)){
    results <- results[results$Estimate != 0,]
    return(list(pen.vars = results[!is.na(results$mfdr),],unpen.vars = unpen.res ))
  }
  
  ### If cutoff is null return 1:number, else return the minimum
  if(is.null(cutoff)) {
    nShow <- min(number, p)
  } else {
    nShow <- min(sum(results$mfdr <= cutoff, na.rm = TRUE), number, p)
  }
  results$Selected <- ifelse(results$Estimate != 0, "*"," ")
  return(list(pen.vars = results[1:nShow,], unpen.vars = unpen.res))
}
