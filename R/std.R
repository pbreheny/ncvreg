std <- function(XX) {
  STD <- .Call("standardize", XX)
  X <- STD[[1]]
  attr(X, "center") <- STD[[2]]
  attr(X, "center") <- STD[[2]]
  attr(X, "scale") <- STD[[3]]
  X
}
