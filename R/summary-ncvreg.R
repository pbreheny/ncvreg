#' Summary method for ncvreg objects
#' 
#' Inferential summaries for `ncvreg` and `ncvsurv` objects based on local marginal false discovery rates.
#' 
#' @param object An `ncvreg` or `ncvsurv` object.
#' @param lambda The regularization parameter value at which inference should be
#'   reported.
#' @param which Alternatively, `lambda` may be specified by index; `which=10`
#'   means: report inference for the 10th value of `lambda` along the
#'   regularization path.  If both `lambda` and `which` are specified, `lambda`
#'   takes precedence.
#' @param number By default, `summary` will provide an inferential summary for
#'   each variable that has been selected (i.e.,  each variable with a nonzero
#'   coefficient). Specifying `number=5`, for example, means that the summary
#'   table will include the 5 features with the lowest mfdr values, regardless
#'   of whether they were selected.  To see all features, `number=Inf`.
#' @param cutoff Alternatively, specifying for example `cutoff=0.3` will report
#'   inference for all features with mfdr under 30%. If both `number` and
#'   `cutoff` are specified, the intersection between both sets of features is
#'   reported.
#' @param sort Should the results be sorted by `mfdr`? (default: TRUE)
#' @param sigma For linear regression models, users can supply an estimate of
#'   the residual standard deviation. The default is to use RSS / DF, where
#'   degrees of freedom are approximated using the number of nonzero
#'   coefficients.
#' @param ... Further arguments; in particular, if you have set `returnX=FALSE`,
#'   you will need to supply `X` and `y` in order to calculate local mFDRs.
#' 
#' @returns An object with S3 class `summary.ncvreg`. The class has its own
#' print method and contains the following list elements:
#' \item{penalty}{The penalty used by `ncvreg` or `ncvsurv`}
#' \item{model}{Either `"linear"`, `"logistic"`, or `"Cox"`.}
#' \item{n}{Number of instances.}
#' \item{p}{Number of regression coefficients (not including the intercept).}
#' \item{lambda}{The `lambda` value at which inference is being reported.}
#' \item{nvars}{The number of nonzero coefficients (again, not including the intercept) at that value of `lambda`.}
#' \item{table}{A table containing estimates, normalized test statistics (z), and an estimate of the local mfdr for each coefficient. The mfdr may be loosely interpreted, in an empirical Bayes sense, as the probability that the given feature is null.}
#' \item{unpen.table}{If there are any unpenalized coefficients, a separate inferential summary is given for them.  Currently, this is based on `lm`/`glm`/`coxph` using the penalized coefficients to provide an offset. This is useful and more or less accurate, but not ideal; we hope to improve the inferential methods for unpenalized variables in the future.}
#'    
#' @author Patrick Breheny <patrick-breheny@@uiowa.edu>
#' 
#' @seealso [ncvreg()], [cv.ncvreg()], [plot.cv.ncvreg()], [local_mfdr()]
#' 
#' @examples
#' # Linear regression --------------------------------------------------
#' data(Prostate)
#' fit <- ncvreg(Prostate$X, Prostate$y)
#' summary(fit, lambda=0.08)
#' 
#' # Logistic regression ------------------------------------------------
#' data(Heart)
#' fit <- ncvreg(Heart$X, Heart$y, family="binomial")
#' summary(fit, lambda=0.05)
#' 
#' # Cox regression -----------------------------------------------------
#' data(Lung)
#' fit <- ncvsurv(Lung$X, Lung$y)
#' summary(fit, lambda=0.1)
#' 
#' # Options ------------------------------------------------------------
#' fit <- ncvreg(Heart$X, Heart$y, family="binomial")
#' summary(fit, lambda=0.08, number=3)
#' summary(fit, lambda=0.08, number=Inf)
#' summary(fit, lambda=0.08, cutoff=0.5)
#' summary(fit, lambda=0.08, number=3, cutoff=0.5)
#' summary(fit, lambda=0.08, number=5, cutoff=0.1)
#' summary(fit, lambda=0.08, number=Inf, sort=FALSE)
#' summary(fit, lambda=0.08, number=3, cutoff=0.5, sort=FALSE)
#' 
#' # If X and y are not returned with the fit, they must be supplied
#' fit <- ncvreg(Heart$X, Heart$y, family="binomial", returnX=FALSE)
#' summary(fit, X=Heart$X, y=Heart$y, lambda=0.08)
#' @export

summary.ncvreg <- function(object, lambda, which, number, cutoff, sort=TRUE, sigma, ...) {
  nvars <- predict(object, type="nvars", lambda=lambda, which=which)
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)", call.=FALSE)
  if (missing(lambda)) lambda <- object$lambda[which]
  custom <- !(missing(number) & missing(cutoff))
  if (inherits(object, 'ncvsurv')) {
    model <- 'Cox'
  } else {
    model <- switch(object$family, gaussian="linear", binomial="logistic", poisson="Poisson")
  }
  local <- local_mfdr(object, lambda, sigma=sigma, ...)
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
  if (!missing(number) && number > nrow(Tab)) number <- nrow(Tab)
  if (missing(number) & missing(cutoff)) {
    Tab <- Tab[Tab$Estimate != 0,]
  } else if (missing(cutoff)) {
    Tab <- Tab[Tab$mfdr <= sort(Tab$mfdr)[number],]
  } else if (missing(number)) {
    Tab <- Tab[Tab$mfdr <= cutoff,]
  } else {
    Tab <- Tab[Tab$mfdr <= min(sort(Tab$mfdr)[number], cutoff),]
  }
  if (sort) Tab <- Tab[order(Tab$mfdr),]
  out$table <- Tab
  out
}

#' @param x        A `summary.ncvreg` object.
#' @param digits   Number of digits past the decimal point to print out. Can be a vector specifying different display digits for
#'   each of the five non-integer printed values.
#' 
#' @rdname summary.ncvreg
#' @export

print.summary.ncvreg <- function(x, digits, ...) {
  digits <- if (missing(digits)) digits <- c(4, 2, 2, 3) else rep(digits, length.out=5)
  caller <- sys.function(sys.parent())
  if (!is.null(caller) && identical(caller, print.summary.cv.ncvreg)) {
    skip_intro <- TRUE
  } else {
    skip_intro <- FALSE
  }
  if (!skip_intro) {
    cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, ", p=", x$p, "\n", sep="")
    cat("At lambda=", formatC(x$lambda, digits[1], format="f"), ":\n", sep="")
    cat("-------------------------------------------------\n")
    if (x$custom) {
      cat("  Features satisfying criteria       : ", nrow(x$table), "\n", sep="")
      cat("  Average mfdr among chosen features : ", formatC(mean(x$table$mfdr), digits=3), "\n", sep="")
    } else {
      cat("  Nonzero coefficients         : ", formatC(x$nvars, width=3), "\n", sep="")
    }
  }
  if (!x$custom) {
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
