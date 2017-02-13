set.seed(1)

#################################################
.test = "ncvreg works for logistic regression" ##
#################################################
X <- matrix(rnorm(500),ncol=10)
b <- rnorm(10)
y <- rnorm(X%*%b) > 0
beta <- glm(y~X,family="binomial")$coef
fit <- ncvreg(X,y,lambda=0,family="binomial",penalty="SCAD",eps=.0001)
scad <- coef(ncvreg(X,y,lambda=0,family="binomial",penalty="SCAD",eps=.0001))
mcp <- coef(ncvreg(X,y,lambda=0, family="binomial",penalty="MCP", eps=.0001))
check(scad, beta,tolerance=.01,check.attributes=FALSE)
check(mcp, beta,tolerance=.01,check.attributes=FALSE)

##############################################
.test = "ncvreg reproduces lasso: binomial" ##
##############################################
require(glmnet)
nlasso <- coef(fit <- ncvreg(X, y, family="binomial", penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, y, family="binomial", lambda=fit$lambda)))
plot(fit, "lambda")
check(nlasso,  glasso, tolerance=.01, check.attributes=FALSE)

##########################################
.test = "logLik() is correct: binomial" ##
##########################################
fit.mle <- glm(y~X, family="binomial")
fit <- ncvreg(X, y, lambda.min=0, family="binomial")
check(logLik(fit)[100],  logLik(fit.mle)[1], check.attributes=FALSE, tol= .001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol= .001)

###############################################
.test = "ncvreg dependencies work: binomial" ##
###############################################

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

##################################################
.test = "cv.ncvreg() options work for binomial" ##
##################################################
n <- 200
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=1) > 0.5

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y, family="binomial")
plot(cvfit, type="all")
print(summary(cvfit))
print(predict(cvfit, type="coefficients"))
print(predict(cvfit, type="vars"))
print(predict(cvfit, type="nvars"))
print(head(predict(cvfit, X=X, "link")))
print(head(predict(cvfit, X=X, "response")))
print(head(predict(cvfit, X=X, "class")))

b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=5) > 0.5
cvfit <- cv.ncvreg(X, y, family="binomial")
plot(cvfit, type="all")

b <- rep(0, 10)
y <- rnorm(n, mean=X%*%b, sd=5) > 0.5
cvfit <- cv.ncvreg(X, y, family="binomial")
plot(cvfit, type="all")
