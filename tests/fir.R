set.seed(1)

#### Permutation ####

X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)

pmfit <- perm.ncvreg(X, y)
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)

pmfit <- perm.ncvreg(X, y>0, family='binomial')
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)

pmfit <- perm.ncvreg(X, rank(y), family='poisson')
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)

#### Permutation of residuals ####
pmfit <- perm.ncvreg(X, y, permute='residuals')
plot(pmfit)

#### Analytic ####
fit <- ncvreg(X, y)
fir(fit)
plot(fir(fit))
plot(pmfit, type="EF")
plot(fir(fit), log.l=TRUE)
