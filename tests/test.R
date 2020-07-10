if (requireNamespace("tinytest", quietly=TRUE)) {
  if (length(unclass(packageVersion("ncvreg"))[[1]]) == 4 | Sys.getenv('R_FORCE_TEST') == 'TRUE') {
    tinytest::test_package("ncvreg", pattern="^[^_].*\\.[rR]$")
  }
}
