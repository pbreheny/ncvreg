if (requireNamespace("tinytest", quietly=TRUE)) {
  if (length(unclass(packageVersion("ncvreg"))[[1]]) == 4) {
    tinytest::test_package("ncvreg", pattern="^[^_].*\\.[rR]$")
  }
}
