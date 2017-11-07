std <- function(X) {
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  STD <- .Call("standardize", X)
  dimnames(STD[[1]]) <- dimnames(X)
  ns <- which(STD[[3]] > 1e-6)
  if (length(ns) == ncol(X)) {
    val <- STD[[1]]
  } else {
    val <- STD[[1]][, ns, drop=FALSE]
  }
  attr(val, "center") <- STD[[2]]
  attr(val, "scale") <- STD[[3]]
  attr(val, "nonsingular") <- ns
  val
}
