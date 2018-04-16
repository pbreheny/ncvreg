summary.ncvsurv <- function(object, lambda, which, ...) {
  nvars <- predict(object, type="nvars", lambda=lambda, which=which)
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)")
  if (missing(lambda)) lambda <- object$lambda[which]
  model <- "Cox"
  local <- local_mfdr(object, lambda, ...)
  val <- list(penalty=object$penalty, model=model, n=object$n, p=nrow(object$beta), lambda=lambda, nvars=nvars, table=local$pen.vars, unpenTable = local$unpen.vars)
  if ("X" %in% names(object)) {
    mFDR <- mfdr(object)
    f <- approxfun(object$lambda, mFDR$EF)
    val$EF = f(lambda)
  }
  structure(val, class="summary.ncvreg")
}
