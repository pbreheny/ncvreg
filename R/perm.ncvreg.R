perm.ncvreg <- function(X, y, ..., N=10, seed, trace=FALSE) {
  if (!missing(seed)) set.seed(seed)
  fit <- ncvreg(X=X, y=y, returnX=TRUE, ...)
  S <- predict(fit, type="nvars")
  p <- ncol(fit$X)
  nlambda <- ncol(fit$beta)

  S.perm <- L.perm <- matrix(NA, N, length(fit$lam))
  for (i in 1:N) {
    if (trace) cat("Starting permutation fit #",i,sep="","\n")

    ## Fit
    if (fit$family=="gaussian") {
      res <- .Call("cdfit_gaussian", fit$X, sample(fit$y), fit$penalty, fit$lambda, 0.001, as.integer(1000), as.double(fit$gamma), fit$penalty.factor, fit$alpha, as.integer(max(S)), as.integer(TRUE))
      b <- matrix(res[[1]], p, nlambda)
      L.perm[i,] <- res[[2]]
    } else if (fit$family=="binomial") {
      res <- .Call("cdfit_binomial", fit$X, sample(fit$y), fit$penalty, fit$lambda, 0.001, as.integer(1000), as.double(fit$gamma), fit$penalty.factor, fit$alpha, as.integer(max(S)), as.integer(TRUE), as.integer(FALSE))
      b <- matrix(res[[2]], p, nlambda)
      L.perm[i,] <- res[[3]]
    }
    
    S.perm[i,] <- apply(b!=0, 2, sum)
  }

  EF <- pmin(apply(S.perm, 2, mean), S)

  FIR <- EF/S
  FIR[S==0] <- 0
  fit <- structure(fit[1:11], class="ncvreg")
  structure(list(EF=EF, S=S, FIR=FIR, fit=fit, loss=apply(L.perm, 2, mean)), class=c("perm.ncvreg", "fir"))
}
