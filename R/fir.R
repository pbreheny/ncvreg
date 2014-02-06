fir <- function(fit) {
  if (fit$family == "binomial") stop("FIR not yet implemented for logistic regression")
  tau <- sqrt(fit$loss)/fit$n
  p <- dim(fit$beta)[1]-1
  S <- predict(fit, type="nvars")
  l <- fit$lam*fit$alpha
  EF <- pmin(p*2*pnorm(-l/tau), S)
  FIR <- EF/S
  FIR[S==0] <- 0
  structure(list(EF=EF, S=S, FIR=FIR), class="fir")
}
