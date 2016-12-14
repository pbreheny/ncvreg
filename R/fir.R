fir <- function(fit) {
  stopifnot(class(fit)[1] %in% c("ncvreg", "ncvsurv"))
  S <- predict(fit, type="nvars")
  if (class(fit)[1] == "ncvsurv") {
    p <- dim(fit$beta)[1]
    R <- apply(fit$W, 2, function(x) rev(cumsum(rev(x))))
    pD <- fit$W/R
    W <- apply(fit$fail*pD*(1-pD), 2, cumsum)
    tau <- sqrt(apply(W, 2, mean))
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
