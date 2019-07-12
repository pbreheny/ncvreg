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
