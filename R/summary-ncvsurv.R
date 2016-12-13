summary.ncvsurv <- function(object, lambda, which, ...) {
  nvars <- predict(object, type="nvars", lambda=lambda, which=which)
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)")
  FIR <- fir(object)
  if (missing(lambda)) lambda <- object$lambda[which]
  f <- approxfun(object$lambda, FIR$EF)
  model <- "Cox"
  val <- list(penalty=object$penalty, model=model, n=object$n, p=nrow(object$beta), lambda=lambda, nvars=nvars, EF=f(lambda))
  structure(val, class="summary.ncvreg")
}
