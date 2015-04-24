source("~/dev/.ncvreg.setup.R")

## Gen data
n <- 40
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- rnorm(n, mean=X%*%b, sd=2)
## yb <- y > .5
## yp <- rpois(n, exp(X%*%b/3))

## Permuting y
pmfit <- perm.ncvreg(X, y)

par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)
plot(pmfit$fit)

## permres
fit <- ncvreg(X, y, returnX=TRUE)
permres.ncvreg(fit, lam=0.2)

## Permuting residuals
pmfit <- perm.ncvreg(X, y, permute="residuals", N=25)
par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)
plot(pmfit$fit)
