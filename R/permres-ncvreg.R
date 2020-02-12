permres.ncvreg <- function(fit, lambda, N=10, seed, trace=FALSE, ...) {
  if (!missing(seed)) set.seed(seed)
  if (!inherits(fit, "ncvreg")) stop("fit is not an 'ncvreg' object", call.=FALSE)
  if (!is.numeric(lambda)) stop("lambda must be numeric", call.=FALSE)
  if (length(lambda) != 1) stop("lambda must be a single number, not a vector", call.=FALSE)
  if (is.null(fit$X)) stop("must run ncvreg with returnX=TRUE", call.=FALSE)

  S <- predict(fit, type="nvars", lambda=lambda)
  y <- fit$y - predict(fit, fit$X, lambda=lambda)
  y <- y - mean(y)
  lam <- c(fit$lambda[fit$lambda > lambda], lambda)
  pfit <- fit.perm.ncvreg(fit, y, lam, N, S, trace)
  S.perm <- pfit$S.perm[, ncol(pfit$S.perm)]
  L.perm <- pfit$L.perm[, ncol(pfit$S.perm)]

  EF <- mean(S.perm, na.rm=TRUE)
  mFDR <- EF/S
  mFDR[S==0] <- 0
  list(EF=EF, S=S, mFDR=mFDR, loss=mean(L.perm, na.rm=TRUE))
}
permres <- function(fit, ...) UseMethod("permres")
