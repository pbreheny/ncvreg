if (interactive()) library(tinytest)

## Works
data(Prostate)
X <- Prostate$X
y <- Prostate$y
boot <- boot_ncvreg(X, y)
# plot(boot)
