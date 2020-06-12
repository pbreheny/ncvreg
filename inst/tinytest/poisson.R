################################################
# ncvreg works for Poisson regression
################################################
n <- 200
p <- 50
X <- matrix(rnorm(n*p), ncol=p)
y <- rpois(n, 1)
beta <- glm(y~X, family="poisson")$coef
scad <- coef(ncvreg(X,y,lambda=0,family="poisson",penalty="SCAD",eps=.0001))
mcp <- coef(ncvreg(X,y,lambda=0, family="poisson",penalty="MCP", eps=.0001))
expect_equivalent(scad, beta,tolerance=.01)
expect_equivalent(mcp, beta,tolerance=.01)

##############################################
# ncvreg reproduces lasso: poisson
##############################################
require(glmnet)
nlasso <- coef(fit <- ncvreg(X, y, family="poisson", penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, y, family="poisson", lambda=fit$lambda)))
plot(fit, "lambda")
expect_equivalent(nlasso,  glasso, tolerance=.01)

################################
# logLik() is correct
################################
fit.mle <- glm(y~X, family="poisson")
fit <- ncvreg(X, y, lambda.min=0, family="poisson")
expect_equivalent(logLik(fit)[100],  logLik(fit.mle)[1], tol= .001)
expect_equivalent(AIC(logLik(fit))[100],  AIC(fit.mle), tol= .001)

##############################################
# ncvreg dependencies work: poisson
##############################################

# Predict
predict(fit, X, 'link')[1:5, 1:5]
predict(fit, X, 'response')[1:5, 1:5]
predict(fit, X, 'coef')[1:5, 1:5]
head(predict(fit, X, 'vars'))
head(predict(fit, X, 'nvars'))

#################################################
# cv.ncvreg() options work for poisson
#################################################
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-1, 1, rep(0, p-2))
y <- rpois(n, exp(X%*%b))

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y, family="poisson")
plot(cvfit, type="all")
summary(cvfit)
head(predict(cvfit, type="coefficients"))
predict(cvfit, type="vars")
predict(cvfit, type="nvars")
head(predict(cvfit, X=X, "link"))
head(predict(cvfit, X=X, "response"))

y <- rpois(n, 1)
cvfit <- cv.ncvreg(X, y, family="poisson")
par(mfrow=c(2,2))
plot(cvfit, type="all")

