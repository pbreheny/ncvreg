summary.ncvreg <- function(object, lambda, which, ...) {
  nvars <- predict(object, type="nvars", lambda=lambda, which=which)
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)")
  FIR <- fir(object)
  if (missing(lambda)) lambda <- object$lambda[which]
  f <- approxfun(object$lambda, FIR$EF)
  model <- switch(object$family, gaussian="linear", binomial="logistic", poisson="Poisson")
  val <- list(penalty=object$penalty, model=model, n=object$n, p=nrow(object$beta)-1, lambda=lambda, nvars=nvars, EF=f(lambda))
  structure(val, class="summary.ncvreg")
}
print.summary.ncvreg <- function(x, digits, ...) {
  digits <- if (missing(digits)) digits <- c(4, 2, 2, 3) else rep(digits, length.out=5)
  cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, ", p=", x$p, "\n", sep="")
  cat("At lambda=", formatC(x$lambda, digits[1], format="f"), ":\n", sep="")
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", x$nvars, "\n", sep="")
  cat("  Expected nonzero coefficients: ", formatC(x$EF, digits=digits[2], format="f"), "\n", sep="")
  cat("  FIR: ", formatC(x$EF/x$nvars, digits=3, format="f"), "\n", sep="")
}
