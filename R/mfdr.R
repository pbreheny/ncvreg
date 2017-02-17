mfdr <- function(fit, X) {
  # Setup
  stopifnot(class(fit)[1] %in% c("ncvreg", "ncvsurv"))
  linear <- class(fit)[1] == "ncvreg" & fit$family == "gaussian"
  S0 <- sum(fit$penalty.factor==0)
  S <- predict(fit, type="nvars") - S0
  if ((class(fit)[1] == "ncvsurv") || (class(fit)[1] == "ncvreg" & fit$family == "binomial")) {
    if (!("X" %in% names(fit))) {
      if (missing(X)) {
        stop("For Cox/GLM models, you must either supply X or run ncvsurv with
returnX=TRUE to calculate an mFDR")
      } else {
        if (class(fit)[1] == "ncvsurv") {
          fit$X <- std(X)[fit$order,]
        } else {
          fit$X <- std(X)
        }
      }
    }
  }

  # Call C functions
  if (class(fit)[1] == "ncvsurv") {
    EF <- .Call("mfdr_cox", fit)
  } else {
    if (fit$family == "binomial") {
      EF <- .Call("mfdr_binomial", fit)
    } else if (fit$family == "gaussian") {
      EF <- .Call("mfdr_gaussian", fit)
    }
  }

  # Calculate rate, return
  EF <- pmin(EF - S0, S)
  mFDR <- EF/S
  mFDR[S==0] <- 0
  df <- data.frame(EF=EF, S=S, mFDR=mFDR)
  rownames(df) <- lamNames(fit$lambda)
  structure(df, class=c("mfdr", "data.frame"))
}
fir <- function(fit, ...) {
  warning("
fir has been deprecated and renamed mfdr; please use mfdr() in the future,\n
as support for fir() is likely to be discontinued at some point.")
  mfdr(fit, ...)
}
