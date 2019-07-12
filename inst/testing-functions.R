# These functions are only used in package testing

check <- function(x, y, check.attributes=FALSE, ...) {
  if (missing(y)) {
    xname <- gsub("()", "", match.call()[2])
    if (x==TRUE) return(TRUE)
    message <- paste0("Problem in ", .test, "\n", xname, " FALSE")
  }
  checkResult <- all.equal(x, y, check.attributes=check.attributes, ...)
  if (class(checkResult)[1]=="logical") return(TRUE)
  xname <- gsub("()", "", match.call()[2])
  yname <- gsub("()", "", match.call()[3])
  if (!exists('.test')) .test <- ''
  message <- paste0("Problem in ", .test, "\n", xname, " not equal to ", yname, "\n", checkResult)
  stop(message, call.=FALSE)
}

median.survfit <- function(x, ...) {
  start.time <- if (!is.null(x$start.time)) x$start.time else min(0, x$time)
  pfun <- function(nused, time, surv, n.risk, n.event, lower, upper, start.time, end.time) {
    minmin <- function(y, x) {
      tolerance <- .Machine$double.eps^0.5
      keep <- (!is.na(y) & y < (0.5 + tolerance))
      if (!any(keep))
        NA
      else {
        x <- x[keep]
        y <- y[keep]
        if (abs(y[1] - 0.5) < tolerance && any(y < y[1]))
          (x[1] + x[min(which(y < y[1]))])/2
        else x[1]
      }
    }
    if (!is.na(end.time)) {
      hh <- ifelse((n.risk - n.event) == 0, 0, n.event/(n.risk * (n.risk - n.event)))
      keep <- which(time <= end.time)
      if (length(keep) == 0) {
        temptime <- end.time
        tempsurv <- 1
        hh <- 0
      } else {
        temptime <- c(time[keep], end.time)
        tempsurv <- c(surv[keep], surv[max(keep)])
        hh <- c(hh[keep], 0)
      }
      n <- length(temptime)
      delta <- diff(c(start.time, temptime))
      rectangles <- delta * c(1, tempsurv[-n])
      varmean <- sum(cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1])
      mean <- sum(rectangles) + start.time
    } else {
      mean <- 0
      varmean <- 0
    }
    med <- minmin(surv, time)
    if (!is.null(upper)) {
      upper <- minmin(upper, time)
      lower <- minmin(lower, time)
      c(nused, max(n.risk), n.risk[1], sum(n.event), sum(mean),
        sqrt(varmean), med, lower, upper)
    } else c(nused, max(n.risk), n.risk[1], sum(n.event), sum(mean), sqrt(varmean), med, 0, 0)
  }
  stime <- x$time
  surv <- x$surv
  plab <- c("records", "n.max", "n.start", "events", "*rmean", "*se(rmean)", "median", paste(x$conf.int, c("LCL", "UCL"), sep = ""))
  ncols <- 9
  if (is.matrix(surv) && !is.matrix(x$n.event)) x$n.event <- matrix(rep(x$n.event, ncol(surv)), ncol = ncol(surv))
  if (is.null(x$strata)) {
    end.time <- max(stime)
    if (is.matrix(surv)) {
      out <- matrix(0, ncol(surv), ncols)
      for (i in 1:ncol(surv)) {
        if (is.null(x$conf.int))
          out[i, ] <- pfun(x$n, stime, surv[, i], x$n.risk,
                           x$n.event[, i], NULL, NULL, start.time, end.time)
        else out[i, ] <- pfun(x$n, stime, surv[, i],
                              x$n.risk, x$n.event[, i], x$lower[, i], x$upper[,
                                                                              i], start.time, end.time)
      }
      dimnames(out) <- list(dimnames(surv)[[2]], plab)
    } else {
      out <- matrix(pfun(x$n, stime, surv, x$n.risk, x$n.event, x$lower, x$upper, start.time, end.time), nrow = 1)
      dimnames(out) <- list(NULL, plab)
    }
  } else {
    nstrat <- length(x$strata)
    stemp <- rep(1:nstrat, x$strata)
    end.time <- last.time <- (rev(stime))[match(1:nstrat, rev(stemp))]
    if (is.matrix(surv)) {
      ns <- ncol(surv)
      out <- matrix(0, nstrat * ns, ncols)
      if (is.null(dimnames(surv)[[2]]))
        dimnames(out) <- list(rep(names(x$strata), ns), plab)
      else {
        cname <- outer(names(x$strata), dimnames(surv)[[2]], paste, sep = ", ")
        dimnames(out) <- list(c(cname), plab)
      }
      k <- 0
      for (j in 1:ns) {
        for (i in 1:nstrat) {
          who <- (stemp == i)
          k <- k + 1
          if (is.null(x$lower))
            out[k, ] <- pfun(x$n[i], stime[who], surv[who, j], x$n.risk[who], x$n.event[who, j], NULL, NULL, start.time, end.time[i])
          else out[k, ] <- pfun(x$n[i], stime[who], surv[who, j], x$n.risk[who], x$n.event[who, j], x$lower[who, j], x$upper[who, j], start.time, end.time[i])
        }
      }
    }
    else {
      out <- matrix(0, nstrat, ncols)
      dimnames(out) <- list(names(x$strata), plab)
      for (i in 1:nstrat) {
        who <- (stemp == i)
        if (is.null(x$lower))
          out[i, ] <- pfun(x$n[i], stime[who], surv[who], x$n.risk[who], x$n.event[who], NULL, NULL, start.time, end.time[i])
        else out[i, ] <- pfun(x$n[i], stime[who], surv[who], x$n.risk[who], x$n.event[who], x$lower[who], x$upper[who], start.time, end.time[i])
      }
    }
  }
  drop(out[,'median'])
}
