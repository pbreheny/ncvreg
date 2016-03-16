set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001, check.attributes=FALSE)}

#### Linear regression ####

# Works
X <- matrix(rnorm(500), 50, 10)
y <- rnorm(50)
fit <- ncvreg(X, y, lambda.min=0, penalty="MCP")

# Equals MLE when lam=0
beta <- lm(y~X)$coef
stopifnot(equal(coef(fit)[,100], beta))

# Different penalties
fit <- ncvreg(X, y, lambda.min=0, penalty="SCAD")
stopifnot(equal(coef(fit)[,100], beta))
fit <- ncvreg(X, y, lambda.min=0, penalty="lasso")
stopifnot(equal(coef(fit)[,100], beta))

# Integers
X <- matrix(rpois(500, 1), 50, 10)
y <- rpois(50, 1)
fit <- ncvreg(X, y)

# Data frame
fit <- ncvreg(as.data.frame(X), y)

# Penalty factor
fit <- ncvreg(as.data.frame(X), y, penalty.factor=c(0:9))

# User lambdas
fit <- ncvreg(as.data.frame(X), y, lambda=c(1, 0.1, 0.01))

# ReturnX
fit <- ncvreg(as.data.frame(X), y, returnX=TRUE)

#### Logistic regression ####

# Works
y <- rbinom(50, 1, 0.5)
fit <- ncvreg(X, y, lambda.min=0, family='binomial')

# Equals MLE when lam=0
beta <- glm(y~X, family='binomial')$coef
stopifnot(equal(coef(fit)[,100], beta))

# y logical
fit <- ncvreg(X, y==1, lambda.min=0, family='binomial')

#### Poisson regression ####

fit <- ncvreg(X, y, lambda.min=0, family='poisson')
