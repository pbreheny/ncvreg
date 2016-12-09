fir <- function(fit) {
  stopifnot(class(fit)[1] %in% c("ncvreg", "ncvsurv"))
  S <- predict(fit, type="nvars")
  if (class(fit)[1] == "ncvsurv") {
    p <- dim(fit$beta)[1]
    R <- apply(fit$W, 2, function(x) rev(cumsum(rev(x))))
    pD <- fit$W/R
    W <- apply(fit$fail*pD*(1-pD), 2, cumsum)
    tau <- apply(W, 2, mean)
  } else {
    p <- dim(fit$beta)[1]-1
    if (fit$family=="gaussian") {
      tau <- sqrt(fit$loss/(fit$n - p + 1))
    } else if (fit$family=="binomial") {
      tau <- sqrt(fit$wMean)
    } else {
      stop(paste0("FIR not yet implemented for family=", fit$family))
    }
  }
  if (all(fit$penalty.factor==1)) {
    l <- fit$lam*fit$alpha
    EF <- pmin(2*p*pnorm(-sqrt(fit$n)*l/tau), S)
  } else {
    S0 <- sum(fit$penalty.factor==0)
    S <- S - S0
    L <- outer(fit$lambda, fit$alpha*fit$penalty.factor)
    EF <- pmin(apply(2*pnorm(-sqrt(fit$n)*L/tau), 1, sum) - S0, S)
  }
  FIR <- EF/S
  FIR[S==0] <- 0
  df <- data.frame(EF=EF, S=S, FIR=FIR)
  rownames(df) <- lamNames(fit$lambda)
  structure(df, class=c("fir", "data.frame"))
}

fir.ncvsurv <- function(fit) {
  X <- fit$X
  d <- fit$fail
  S <- predict(fit, type="nvars")
  l <- fit$lambda*fit$alpha
  pf <- which(fit$penalty.factor != 0)
  p <- 40
  t <- 4
  EF <- WW <- NULL
  for (i in 1:length(l)){
    ### Calculate W matrix
    R <- rev(cumsum(rev(fit$W[,i])))
    Rk <- 1/R
    Rk2 <- 1/(R^2)
    RR <- cumsum(Rk)
    RR2 <- cumsum(Rk2)
    W <- diag(fit$W[,i]*RR - (fit$W[,i]^2)*RR2)

    ### Calculate x'Wx
    for (j in 1:ncol(fit$X)){
      WW[j] <- t(X[,j]) %*% W %*% X[,j]
    }

    ### Calculate expected false inclusions
    EF[i] <- 2*sum(pnorm(-(nrow(X)*l[i])/sqrt(WW)))
  }
  EF <- pmin(EF, S)
  FIR <- EF/S
  FIR[S==0] <- 0
  df <- data.frame(EF=EF, S=S, FIR=FIR)
  structure(df, class=c("fir", "data.frame"))
}
