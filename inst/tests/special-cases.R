############################################
.test = "ncvreg handles constant columns" ##
############################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
X[, 3:5] <- 0
y <- rnorm(n)
fit <- ncvreg(X, y)
y <- runif(n) > .5
fit <- ncvreg(X, y, family="binomial")
y <- rpois(n, 1)
fit <- ncvreg(X, y, family="poisson")

############################################
.test = "ncvsurv handles constant columns" ##
############################################
y <- cbind(rexp(n), rbinom(n, 1, prob=0.8))
fit <- ncvsurv(X, y)

#########################################
.test = "ncvreg handles single lambda" ##
#########################################

n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)

y <- rnorm(n)
fit <- ncvreg(X, y, lambda=0.1)
coef(fit)

y <- rbinom(n, 1, 0.5)
fit <- ncvreg(X, y, lambda=0.02, family="binomial")
coef(fit)

y <- rpois(n, 1)
fit <- ncvreg(X, y, lambda=0.05, family="poisson")
coef(fit)

##########################################
.test = "ncvsurv handles single lambda" ##
##########################################

y <- cbind(rexp(n, 1), 1)
fit <- ncvsurv(X, y, lambda=0.05)
coef(fit)
