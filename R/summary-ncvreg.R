summary.ncvreg <- function(object, lambda, which, number, cutoff, ...) {
  nvars <- predict(object, type="nvars", lambda=lambda, which=which)
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)")
  if (missing(lambda)) lambda <- object$lambda[which]
  if (missing(number)) number <- NULL
  if (missing(cutoff)) cutoff <- NULL
  custom <- !(is.null(number) & is.null(cutoff))
  if ('ncvsurv' %in% class(object)) {
    model <- 'Cox'
  } else {
    model <- switch(object$family, gaussian="linear", binomial="logistic", poisson="Poisson")
  }
  local <- local_mfdr(object, lambda, number, cutoff, ...)
  structure(
    list(penalty=object$penalty, model=model, n=object$n, p=nrow(object$beta)-1, lambda=lambda, nvars=nvars, table=local$pen.vars, unpenTable = local$unpen.vars, custom=custom),
    class='summary.ncvreg'
  )
}
print.summary.ncvreg <- function(x, digits, ...) {
  digits <- if (missing(digits)) digits <- c(4, 2, 2, 3) else rep(digits, length.out=5)
  cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, ", p=", x$p, "\n", sep="")
  cat("At lambda=", formatC(x$lambda, digits[1], format="f"), ":\n", sep="")
  cat("-------------------------------------------------\n")
  if (x$custom) {
    cat("  Features satisfying criteria       : ", nrow(x$table), "\n", sep="")
    cat("  Average mfdr among chosen features : ", formatC(mean(x$table$mfdr), digits=3), "\n", sep="")
  } else {
    cat("  Nonzero coefficients         : ", x$nvars, "\n", sep="")
    if (x$nvars > 0) {
      cat("  Expected nonzero coefficients: ", formatC(sum(x$table$mfdr), digits=3), '\n', sep="")
      space <- substr('                                              ', 1, 5-nchar(nrow(x$table)))
      cat("  Average mfdr (", nrow(x$table), " features)", space, ": ", formatC(mean(x$table$mfdr), digits=3), '\n', sep="")
    }
  }
  cat("\n")
  x$table$mfdr <- format.pval(x$table$mfdr, eps=1e-4)
  print(x$table, digits=digits)
  if (!is.null(x$unpenTable)) {
    x$unpenTable$p.value <- format.pval(x$unpenTable$p.value, eps=1e-4)
    cat("\nUnpenalized variables:\n")
    print(x$unpenTable, digits=digits)
  }
}
