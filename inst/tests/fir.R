require(ncvreg)
set.seed(1)

#### Permutation ####

X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)

pmfit <- perm.ncvreg(X, y)
pmfit <- perm.ncvreg(X, y>0, family='binomial')
pmfit <- perm.ncvreg(X, rank(y), family='poisson')

#### Permutation of residuals ####
pmfit <- perm.ncvreg(X, y, permute='residuals')

#### Analytic ####
fit <- ncvreg(X, y)
fir(fit)
