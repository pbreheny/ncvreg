n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)

.test <- "gaussian: single lambda"
y <- rnorm(n)
fit <- ncvreg(X, y, lambda=0.1)
coef(fit)

.test <- "binomial: single lambda"
y <- rbinom(n, 1, 0.5)
fit <- ncvreg(X, y, lambda=0.02, family="binomial")
coef(fit)

.test <- "Poisson: single lambda"
y <- rpois(n, 1)
fit <- ncvreg(X, y, lambda=0.05, family="poisson")
coef(fit)

.test <- "Cox: single lambda"
y <- cbind(rexp(n, 1), 1)
fit <- ncvsurv(X, y, lambda=0.05)
coef(fit)
