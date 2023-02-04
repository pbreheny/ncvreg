lam_names <- function(l) {
  if (length(l) > 1) {
    d <- ceiling(-log10(-max(diff(l))))
    d <- min(max(d,4), 10)
  } else {
    d <- 4
  }
  formatC(l, format="f", digits=d)
}
