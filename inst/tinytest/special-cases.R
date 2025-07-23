# --- ncvreg handles constant columns -------------------------------
n <- 50
p <- 10
X <- matrix(rnorm(n * p), ncol = p)
X[, 3:5] <- 0
y <- rnorm(n)
fit <- ncvreg(X, y)
y <- runif(n) > .5
fit <- ncvreg(X, y, family = "binomial")
y <- rpois(n, 1)
fit <- ncvreg(X, y, family = "poisson")

# --- ncvsurv handles constant columns ------------------------------
y <- cbind(rexp(n), rbinom(n, 1, prob = 0.8))
fit <- ncvsurv(X, y)
