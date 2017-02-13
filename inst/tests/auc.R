require(survival)
#source("~/dev/.ncvreg.setup.R")
n <- 500
X <- matrix(rnorm(n*10), n, 10)
y <- Surv(rexp(n, rate=exp(X[,1])), rbinom(n, 1, prob=0.8))
cvfit <- cv.ncvsurv(X, y, lambda.min=0, returnY=TRUE)
a <- AUC(cvfit)
a[length(a)]
fit <- coxph(y~X)
survConcordance(y~predict(fit))$concordance
