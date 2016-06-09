std <- function(X) {
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  STD <- .Call("standardize", X)
  val <- STD[[1]]
  attr(val, "center") <- STD[[2]]
  attr(val, "center") <- STD[[2]]
  attr(val, "scale") <- STD[[3]]
  val
}
