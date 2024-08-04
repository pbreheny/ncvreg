if (interactive()) library(tinytest)

## Works
data(Prostate)
X <- Prostate$X
y <- Prostate$y
boot <- boot_ncvreg(X, y)
confint(boot)

## Specify parm
tinytest::expect_equal(confint(boot, parm = 1:4)$variable, colnames(X)[1:4])
tinytest::expect_equal(confint(boot, parm = "age")$variable, "age")

# Change level
confint(boot, level = 0.9)

## Warning if NA in draws
boot$draws[1,1] <- NA
tinytest::expect_warning(confint(boot))

## Expect error
tinytest::expect_error(confint(boot, level = 1.1))
tinytest::expect_error(confint(boot, level = 2))
tinytest::expect_error(confint(boot, parm = mean))

## Checks
tinytest::expect_equal(class(confint(boot)), "data.frame")
tinytest::expect_equal(colnames(confint(boot)), c("variable", "estimate", "lower", "upper"))

