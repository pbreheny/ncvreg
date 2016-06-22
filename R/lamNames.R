lamNames <- function(l) {
  d <- ceiling(-log10(-max(diff(l))))
  d <- min(max(d,4), 10)
  formatC(l, format="f", digits=d)
}
