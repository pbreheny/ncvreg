summary.ncvreg <- function(object, lambda, which, ...) {
  nvars <- predict(object, type="nvars", lambda=lambda, which=which)
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)")
  if (missing(lambda)) lambda <- object$lambda[which]
  model <- switch(object$family, gaussian="linear", binomial="logistic", poisson="Poisson")
  local <- local_mfdr(object, lambda, ...)
  val <- list(penalty=object$penalty, model=model, n=object$n, p=nrow(object$beta)-1, lambda=lambda, nvars=nvars, table=local$pen.vars, unpenTable = local$unpen.vars)
  if ("X" %in% names(object) || object$family=="gaussian") {
    mFDR <- mfdr(object)
    f <- approxfun(object$lambda, mFDR$EF)
    val$EF = f(lambda)
  }
  structure(val, class="summary.ncvreg")
}
print.summary.ncvreg <- function(x, digits, ...) {
  digits <- if (missing(digits)) digits <- c(4, 2, 2, 3) else rep(digits, length.out=5)
  cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, ", p=", x$p, "\n", sep="")
  cat("At lambda=", formatC(x$lambda, digits[1], format="f"), ":\n", sep="")
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", x$nvars, "\n", sep="")
  if ("EF" %in% names(x)) {
    cat("  Expected nonzero coefficients: ", formatC(x$EF, digits=digits[2], format="f"), "\n", sep="")
    cat("  mFDR: ", formatC(x$EF/x$nvars, digits=3, format="f"), "\n", sep="")
    if (nrow(x$table) == x$nvars & all(x$table$Estimate != 0)) {
      cat("  (local) Expected nonzero coefficients: ", formatC(sum(x$table$mfdr), digits=3), '\n', sep="")
      cat("  (local) mfdr: ", formatC(mean(x$table$mfdr), digits=3), '\n', sep="")
    }
    cat("\n")
  }
  x$table$mfdr <- format.pval(x$table$mfdr, eps=1e-4)
  print(x$table, digits=digits)
  if (!is.null(x$unpenTable)) {
    x$unpenTable$mfdr <- format.pval(x$unpenTable$mfdr, eps=1e-4)
    print(x$unpenTable, digits=digits)
  }
}
