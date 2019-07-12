summary.ncvreg <- function(object, lambda, which, number, cutoff, ...) {
  nvars <- predict(object, type="nvars", lambda=lambda, which=which)
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)")
  if (missing(lambda)) lambda <- object$lambda[which]
  custom <- !(missing(number) & missing(cutoff))
  if ('ncvsurv' %in% class(object)) {
    model <- 'Cox'
  } else {
    model <- switch(object$family, gaussian="linear", binomial="logistic", poisson="Poisson")
  }
  local <- local_mfdr(object, lambda, ...)
  out <- structure(list(penalty=object$penalty, model=model, n=object$n, p=nrow(object$beta)-1, lambda=lambda,
                        custom=custom, nvars=nvars),
                   class='summary.ncvreg')
  if (is.null(local$unpen)) {
    Tab <- local
  } else {
    Tab <- local$pen.vars
    out$unpenTable <- local$unpen.vars
  }
  
  # Sort/subset table
  Tab <- Tab[order(Tab$mfdr),]
  if (missing(number) & missing(cutoff)) {
    Tab <- Tab[Tab$Estimate != 0,]
  } else if (missing(cutoff)) {
    Tab <- Tab[1:min(number, nrow(Tab)),]
  } else if (missing(number)) {
    Tab <- Tab[1:min(sum(Tab$mfdr <= cutoff), nrow(Tab)),]
  } else {
    Tab <- Tab[1:min(sum(Tab$mfdr <= cutoff), number, nrow(Tab)),]
  }
  out$table <- Tab
  out
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
    cat("  Nonzero coefficients         : ", formatC(x$nvars, width=3), "\n", sep="")
    if (x$nvars > 0) {
      cat("  Expected nonzero coefficients: ", formatC(sum(x$table$mfdr), digits=2, format='f', width=6), '\n', sep="")
      space <- substr('                                              ', 1, 5-nchar(nrow(x$table)))
      cat("  Average mfdr (", nrow(x$table), " features)", space, ": ", formatC(mean(x$table$mfdr), digits=3, format='f', width=7), '\n', sep="")
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
