se.ncvsurv <-function (y, eta, B = 100) {
  cve <- matrix(NA, B, ncol(eta))
  for (b in 1:B) {
    eta.b <- apply(eta, 2, sample, replace=TRUE)
    cve[b,] <- loss.ncvsurv(y, eta.b)
  }
  apply(cve, 2, sd)
}
