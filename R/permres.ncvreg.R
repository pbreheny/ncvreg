permres.ncvreg <- function(fit, lambda, N=10, seed, trace=FALSE, ...) {
  if (!missing(seed)) set.seed(seed)
  if (class(fit)!="ncvreg") stop("fit is not an 'ncvreg' object")
  if (!is.numeric(lambda)) stop("lambda must be numeric")
  if (length(lambda) != 1) stop("lambda must be a single number, not a vector")
  if (is.null(fit$X)) stop("must run ncvreg with returnX=TRUE")

  S <- predict(fit, type="nvars", lambda=lambda)
  y <- fit$y - predict(fit, fit$X, lambda=lambda)
  y <- y - mean(y)
  lam <- c(fit$lambda[fit$lambda > lambda], lambda)
  pfit <- fit.perm.ncvreg(fit, y, lam, N, S, trace)
  S.perm <- pfit$S.perm[,ncol(pfit$S.perm)]
  L.perm <- pfit$L.perm[,ncol(pfit$S.perm)]

  EF <- mean(S.perm, na.rm=TRUE)
  FIR <- EF/S
  FIR[S==0] <- 0
  list(EF=EF, S=S, FIR=FIR, loss=mean(L.perm, na.rm=TRUE))
}
permres <- function(fit, ...) UseMethod("permres")
