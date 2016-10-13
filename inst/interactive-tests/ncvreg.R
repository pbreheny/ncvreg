#source("~/dev/.ncvreg.setup.R")
check <- function(x, y, ...) {
  if (missing(y)) {
    xname <- gsub("()", "", match.call()[2])
    if (x==TRUE) return(TRUE)
    message <- paste0("Problem in ", .test, "\n", xname, " FALSE")
  }
  checkResult <- all.equal(x, y, ...)
  if (class(checkResult)[1]=="logical") return(TRUE)
  xname <- gsub("()", "", match.call()[2])
  yname <- gsub("()", "", match.call()[3])
  message <- paste0("Problem in ", .test, "\n", xname, " not equal to ", yname, "\n", checkResult)
  stop(message, call.=FALSE)
}

###############################################
.test = "ncvreg works for linear regression" ##
###############################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- rnorm(p)
y <- rnorm(n, X%*%b)
beta <- lm(y~X)$coef
scad <- coef(ncvreg(X,y,lambda=0,penalty="SCAD",eps=.0001))
mcp <- coef(ncvreg(X,y,lambda=0,penalty="MCP",eps=.0001))
check(scad, beta,tolerance=.01,check.attributes=FALSE)
check(mcp, beta,tolerance=.01,check.attributes=FALSE)

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

################################################
.test = "ncvreg works for Poisson regression" ##
################################################
X <- matrix(rnorm(500),ncol=10)
b <- rnorm(10)
y <- rnorm(X%*%b) > 0
beta <- glm(y~X,family="poisson")$coef
scad <- coef(ncvreg(X,y,lambda=0,family="poisson",penalty="SCAD",eps=.0001))
mcp <- coef(ncvreg(X,y,lambda=0, family="poisson",penalty="MCP", eps=.0001))
check(scad, beta,tolerance=.01,check.attributes=FALSE)
check(mcp, beta,tolerance=.01,check.attributes=FALSE)

####################################
.test = "ncvreg reproduces lasso" ##
####################################
require(glmnet)
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
par(mfrow=c(3,2))

y <- rnorm(n)
nlasso <- coef(fit <- ncvreg(X, y, penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, y, lambda=fit$lambda)))
plot(fit, "lambda")
check(nlasso,  glasso, tolerance=.01, check.attributes=FALSE)

yy <- runif(n) > .5
nlasso <- coef(fit <- ncvreg(X, yy, family="binomial", penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, yy, family="binomial", lambda=fit$lambda)))
plot(fit, "lambda")
check(nlasso,  glasso, tolerance=.01, check.attributes=FALSE)

yy <- rpois(n, 1)
nlasso <- coef(fit <- ncvreg(X, yy, family="poisson", penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, yy, family="poisson", lambda=fit$lambda)))
plot(fit, "lambda")
check(nlasso,  glasso, tolerance=.01, check.attributes=FALSE)

################################
.test = "logLik() is correct" ##
################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)

y <- rnorm(n)
fit.mle <- lm(y~X)
fit <- ncvreg(X, y, lambda.min=0)
check(logLik(fit)[100],  logLik(fit.mle)[1], check.attributes=FALSE, tol= .001)
check(AIC(fit)[100],  AIC(fit.mle), check.attributes=FALSE, tol= .001)

yy <- runif(n) > .5
fit.mle <- glm(yy~X, family="binomial")
fit <- ncvreg(X, yy, lambda.min=0, family="binomial")
check(logLik(fit)[100],  logLik(fit.mle)[1], check.attributes=FALSE, tol= .001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol= .001)

yy <- rpois(n, 1)
fit.mle <- glm(yy~X, family="poisson")
fit <- ncvreg(X, yy, lambda.min=0, family="poisson")
check(logLik(fit)[100],  logLik(fit.mle)[1], check.attributes=FALSE, tol= .001)
check(AIC(logLik(fit))[100],  AIC(fit.mle), check.attributes=FALSE, tol= .001)

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

######################################
.test = "cv.ncvreg() seems to work" ##
######################################
n <- 40
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- rnorm(n, mean=X%*%b, sd=2)
yb <- y > .5
yp <- rpois(n, exp(X%*%b/3))
par(mfrow=c(3,2))

require(glmnet)
gcvfit <- cv.glmnet(X, y, nfolds=n)
plot(gcvfit)
ncvfit <- cv.ncvreg(X, y, penalty="lasso", lambda=gcvfit$lambda, nfolds=n)
plot(ncvfit)
gcvfit <- cv.glmnet(X, yb, family="binomial", nfolds=n)
plot(gcvfit)
ncvfit <- cv.ncvreg(X, yb, family="binomial", penalty="lasso", lambda=gcvfit$lambda, nfolds=n)
plot(ncvfit)
cvfit <- cv.glmnet(X, yp, family="poisson")
plot(cvfit)
cvfit <- cv.ncvreg(X, yp, family="poisson", penalty="lasso", lambda=cvfit$lambda)
plot(cvfit)

#########################################
.test = "penalty.factor seems to work" ##
#########################################

##################################################
.test = "cv.ncvreg() options work for gaussian" ##
##################################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=1)

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y)
plot(cvfit, type="all")
print(summary(cvfit))
print(predict(cvfit, "coefficients"))
print(predict(cvfit, "vars"))
print(predict(cvfit, "nvars"))

b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.ncvreg(X, y)
plot(cvfit, type="all")

b <- rep(0, 10)
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.ncvreg(X, y)
plot(cvfit, type="all")

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
print(predict(cvfit, "coefficients"))
print(predict(cvfit, "vars"))
print(predict(cvfit, "nvars"))
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

#################################################
.test = "cv.ncvreg() options work for poisson" ##
#################################################
n <- 200
p <- 50
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-1, 1, rep(0, p-2))
y <- rpois(n, exp(X%*%b))

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y, family="poisson")
plot(cvfit, type="all")
print(summary(cvfit))
print(predict(cvfit, "coefficients"))
print(predict(cvfit, "vars"))
print(predict(cvfit, "nvars"))
print(head(predict(cvfit, X=X, "link")))
print(head(predict(cvfit, X=X, "response")))

y <- rpois(n, 1)
cvfit <- cv.ncvreg(X, y, family="poisson")
par(mfrow=c(2,2))
plot(cvfit, type="all")

##############################################
.test = "cv.ncvreg() return LP array works" ##
##############################################
n <- 100
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-3, 3, rep(0, 8))

y <- rnorm(n, mean=X%*%b, sd=1)
cvfit <- cv.ncvreg(X, y, returnY=TRUE)
cve <- apply(cvfit$Y - y, 2, crossprod)/n
check(cve, cvfit$cve, check.attributes=FALSE, tol= .001)

y <- rnorm(n, mean=X%*%b) > 0
cvfit <- cv.ncvreg(X, y, family='binomial', returnY=TRUE)
pe <- apply((cvfit$Y>0.5)!=y, 2, mean)
check(pe, cvfit$pe, check.attributes=FALSE, tol= .001)

####################################
.test = "standardize=FALSE works" ##
####################################
# n <- 200
# p <- 10
# X <- matrix(rnorm(n*p), ncol=p)
# b <- c(-3, 3, rep(0, 8))
# y <- rnorm(n, mean=X%*%b, sd=1)
# X[,1] <- 5*X[,1]

# fit1 <- ncvreg(X, y, penalty="lasso", standardize=FALSE)
# fit2 <- glmnet(X, y, lambda=b1$lam, standardize=FALSE, intercept=FALSE)
# check(fit1$beta,  as.matrix(coef(fit2))[-1,], check.attributes=FALSE, tol= .001)
#
# yy <- y > 0
# b1 <- ncvreg(X, yy, "binomial", "lasso", standardize=FALSE)
# b2 <- glmnet(X, yy, "binomial", lambda=b1$lam, standardize=FALSE, intercept=TRUE)
# check(b1$beta,  as.matrix(coef(b2)), check.attributes=FALSE, tol= .001)
#
# b1 <- ncvreg(X, yy, "poisson", "lasso", standardize=FALSE)
# b2 <- glmnet(X, yy, "poisson", lambda=b1$lam, standardize=FALSE, intercept=TRUE)
# check(b1$beta,  as.matrix(coef(b2)), check.attributes=FALSE, tol= .001)
#
# b1 <- ncvreg(X, y, penalty="MCP", standardize=FALSE)
# b2 <- ncvreg(X, yy, "binomial", penalty="MCP", standardize=FALSE)
# b3 <- ncvreg(X, y, penalty="SCAD", standardize=FALSE)
# b4 <- ncvreg(X, yy, "binomial", penalty="SCAD", standardize=FALSE)
#
# cvfit1 <- cv.ncvreg(X, y, penalty="lasso", standardize=FALSE)
# cvfit2 <- cv.ncvreg(X, y, penalty="lasso")
