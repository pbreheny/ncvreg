if (interactive()) library(tinytest)

## Works
data(Prostate)
X <- Prostate$X
y <- Prostate$y
boot <- boot_ncvreg(X, y)
confint(boot)

## Warning if NA in draws
boot$draws[1,1] <- NA
tinytest::expect_warning(confint(boot))

## Checks
tinytest::expect_equal(class(confint(boot)), "data.frame")
tinytest::expect_equal(colnames(confint(boot)), c("variable", "estimate", "lower", "upper"))
