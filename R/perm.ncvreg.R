perm.ncvreg <- function(X, y, ..., permute=c("outcome", "residuals"), N=10, seed, trace=FALSE) {
  permute <- match.arg(permute)
  if (!missing(seed)) set.seed(seed)
  fit <- ncvreg(X=X, y=y, returnX=TRUE, ...)
  if (fit$family != 'gaussian' & permute=='residuals') stop("Cannot permute residuals with family = ", fit$family)
  S <- predict(fit, type="nvars")

  if (permute=="outcome") {
    pfit <- fit.perm.ncvreg(fit, fit$y, fit$lambda, N, max(S), trace)
    S.perm <- pfit$S.perm
    L.perm <- pfit$L.perm
    EF <- pmin(apply(S.perm, 2, mean), S)
    FIR <- EF/S
    FIR[S==0] <- 0
    loss <- apply(L.perm, 2, mean)
  } else {
    n.l <- length(fit$lambda)
    EF <- FIR <- loss <- numeric(n.l)
    for (i in 1:n.l) {
      pres <- permres(fit, fit$lambda[i], N=N, seed=seed, trace=trace)
      EF[i] <- pres$EF
      FIR[i] <- pres$FIR
      loss[i] <- pres$loss
      names(EF) <- names(FIR) <- names(loss) <- fit$lambda
    }
  }

  fit <- structure(fit[1:11], class="ncvreg") ## Don't return X, y, etc.
  structure(list(EF=EF, S=S, FIR=FIR, fit=fit, loss=loss), class=c("perm.ncvreg", "fir"))
}
fit.perm.ncvreg <- function(fit, y, lam, N, maxdf, trace) {
  n.l <- length(lam)
  p <- ncol(fit$X)
  S.perm <- L.perm <- matrix(NA, N, n.l)
  for (i in 1:N) {
    if (trace) cat("Starting permutation fit #", i, sep="","\n")
    if (fit$family=="gaussian") {
      res <- .Call("cdfit_gaussian", fit$X, sample(y), fit$penalty, lam, 0.001, as.integer(1000), as.double(fit$gamma), fit$penalty.factor, fit$alpha, as.integer(maxdf), as.integer(TRUE))
      b <- matrix(res[[1]], p, n.l)
      ind <- is.na(res[[3]])
      L.perm[i,] <- res[[2]]
      L.perm[i,ind] <- NA
    } else if (fit$family=="binomial") {
      res <- .Call("cdfit_binomial", fit$X, sample(y), fit$penalty, lam, 0.001, as.integer(1000), as.double(fit$gamma), fit$penalty.factor, fit$alpha, as.integer(maxdf), as.integer(TRUE), as.integer(FALSE))
      b <- matrix(res[[2]], p, n.l)
      ind <- is.na(res[[4]])
      L.perm[i,] <- res[[3]]
      L.perm[i,ind] <- NA
    } else if (fit$family=="poisson") {
      res <- .Call("cdfit_poisson", fit$X, sample(y), fit$penalty, lam, 0.001, as.integer(1000), as.double(fit$gamma), fit$penalty.factor, fit$alpha, as.integer(maxdf), as.integer(TRUE), as.integer(FALSE))
      b <- matrix(res[[2]], p, n.l)
      ind <- is.na(res[[4]])
      L.perm[i,] <- res[[3]]
      L.perm[i,ind] <- NA
    }
    S.perm[i,] <- apply(b[,drop=FALSE]!=0, 2, sum)
    S.perm[i,ind] <- NA
  }
  list(S.perm=S.perm, L.perm=L.perm)
}
