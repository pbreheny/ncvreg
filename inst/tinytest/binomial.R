suppressPackageStartupMessages(library(glmnet))

# ncvreg works for logistic regression
X <- matrix(rnorm(500),ncol=10)
b <- rnorm(10)
y <- rnorm(X%*%b) > 0
beta <- glm(y~X,family="binomial")$coef
fit <- ncvreg(X,y,lambda=0,family="binomial",penalty="SCAD",eps=.0001)
scad <- coef(ncvreg(X,y,lambda=0,family="binomial",penalty="SCAD",eps=.0001))
mcp <- coef(ncvreg(X,y,lambda=0, family="binomial",penalty="MCP", eps=.0001))
expect_equivalent(scad, beta,tolerance=.01)
expect_equivalent(mcp, beta,tolerance=.01)

# ncvreg reproduces lasso: binomial
nlasso <- coef(fit <- ncvreg(X, y, family="binomial", penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, y, family="binomial", lambda=fit$lambda)))
plot(fit, "lambda")
expect_equivalent(nlasso,  glasso, tolerance=.01)

# logLik() is correct: binomial
fit.mle <- glm(y~X, family="binomial")
fit <- ncvreg(X, y, lambda.min=0, family="binomial")
expect_equivalent(logLik(fit)[100],  logLik(fit.mle)[1], tol= .001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol= .001)

# ncvreg dependencies work: binomial

# Predict
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'class')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')

# y logical
fit <- ncvreg(X, y==1, lambda.min=0, family='binomial')

# Summary
summary(fit, which=10)
summary(fit, lam=0.05)

# cv.ncvreg() options work for binomial
n <- 200
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=1) > 0.5

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y, family="binomial")
plot(cvfit, type="all")
summary(cvfit)
predict(cvfit, type="coefficients")
predict(cvfit, type="vars")
predict(cvfit, type="nvars")
head(predict(cvfit, X=X, "link"))
head(predict(cvfit, X=X, "response"))
head(predict(cvfit, X=X, "class"))

b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=5) > 0.5
cvfit <- cv.ncvreg(X, y, family="binomial")
plot(cvfit, type="all")

b <- rep(0, 10)
y <- rnorm(n, mean=X%*%b, sd=5) > 0.5
cvfit <- cv.ncvreg(X, y, family="binomial")
plot(cvfit, type="all")
