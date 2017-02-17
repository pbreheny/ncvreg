################################################
.test = "ncvreg works for Poisson regression" ##
################################################
n <- 200
p <- 50
X <- matrix(rnorm(n*p), ncol=p)
y <- rpois(n, 1)
beta <- glm(y~X, family="poisson")$coef
scad <- coef(ncvreg(X,y,lambda=0,family="poisson",penalty="SCAD",eps=.0001))
mcp <- coef(ncvreg(X,y,lambda=0, family="poisson",penalty="MCP", eps=.0001))
check(scad, beta,tolerance=.01,check.attributes=FALSE)
check(mcp, beta,tolerance=.01,check.attributes=FALSE)

##############################################
.test = "ncvreg reproduces lasso: poisson" ##
##############################################
require(glmnet)
nlasso <- coef(fit <- ncvreg(X, y, family="poisson", penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, y, family="poisson", lambda=fit$lambda)))
plot(fit, "lambda")
check(nlasso,  glasso, tolerance=.01, check.attributes=FALSE)

################################
.test = "logLik() is correct" ##
################################
fit.mle <- glm(y~X, family="poisson")
fit <- ncvreg(X, y, lambda.min=0, family="poisson")
check(logLik(fit)[100],  logLik(fit.mle)[1], check.attributes=FALSE, tol= .001)
check(AIC(logLik(fit))[100],  AIC(fit.mle), check.attributes=FALSE, tol= .001)

##############################################
.test = "ncvreg dependencies work: poisson" ##
##############################################

# Predict
predict(fit, X, 'link')[1:5, 1:5]
predict(fit, X, 'response')[1:5, 1:5]
predict(fit, X, 'coef')[1:5, 1:5]
head(predict(fit, X, 'vars'))
head(predict(fit, X, 'nvars'))

#################################################
.test = "cv.ncvreg() options work for poisson" ##
#################################################
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-1, 1, rep(0, p-2))
y <- rpois(n, exp(X%*%b))

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y, family="poisson")
plot(cvfit, type="all")
print(summary(cvfit))
head(predict(cvfit, type="coefficients"))
print(predict(cvfit, type="vars"))
print(predict(cvfit, type="nvars"))
print(head(predict(cvfit, X=X, "link")))
print(head(predict(cvfit, X=X, "response")))

y <- rpois(n, 1)
cvfit <- cv.ncvreg(X, y, family="poisson")
par(mfrow=c(2,2))
plot(cvfit, type="all")

