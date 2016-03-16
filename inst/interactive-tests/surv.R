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
require(survival)
op <- par(no.readonly = TRUE)

##################################################################
.test = "ncvsurv works for simple cox regression, no censoring" ##
##################################################################
n <- 10
X <- matrix(runif(n), n, 1)
y <- Surv(sort(rexp(n, X), decreasing=TRUE), rep(1,n))
X <- X - mean(X)
beta <- coef(coxph(y ~ X))
eta <- X%*%beta
w <- exp(eta)/cumsum(exp(eta))
sum(X * (1-w))
scad <- coef(ncvsurv(X, y, lambda=0, penalty="SCAD", eps=.0001))
mcp <- coef(ncvsurv(X, y, lambda=0, penalty="MCP", eps=.0001))
check(scad, beta, tolerance=.01, check.attributes=FALSE)
check(mcp, beta, tolerance=.01,check.attributes=FALSE)

####################################################
.test = "ncvsurv works for simple cox regression" ##
####################################################
n <- 30
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
y[which.max(y[,1]), 2] <- 0
y[which.min(y[,1]), 2] <- 0
X <- matrix(rnorm(n), n, 1)
beta <- coef(coxph(y ~ X))
scad <- coef(ncvsurv(X, y, lambda=0, penalty="SCAD", eps=.0001))
mcp <- coef(ncvsurv(X, y, lambda=0,penalty="MCP",eps=.0001))
check(scad, beta, tolerance=.01, check.attributes=FALSE)
check(mcp,beta, tolerance=.01, check.attributes=FALSE)

#############################################
.test = "ncvsurv works for cox regression" ##
#############################################
n <- 100
p <- 10
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
beta <- coef(coxph(y ~ X))
scad <- coef(ncvsurv(X, y, lambda=0, penalty="SCAD", eps=.0001))
mcp <- coef(ncvsurv(X, y, lambda=0,penalty="MCP",eps=.0001))
check(scad, beta, tolerance=.01, check.attributes=FALSE)
check(mcp, beta, tolerance=.01, check.attributes=FALSE)

#######################################
.test = "ncvsurv agrees with coxnet" ##
#######################################
require(glmnet)
n <- 100
p <- 25
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
par(mfrow=c(2,1))

nlasso <- coef(nfit <- ncvsurv(X, y, penalty="lasso", lambda.min=0.001))
plot(nfit, log=TRUE)
glasso <- as.matrix(coef(gfit <- glmnet(X, y, family="cox", lambda=nfit$lambda)))
plot(gfit, "lambda")
check(nlasso, glasso, tolerance=.02, check.attributes=FALSE)

check(predict(nfit, X, "link"), predict(gfit, X, type="link"), tolerance=.01, check.attributes=FALSE)
check(predict(nfit, X, "response"), predict(gfit, X, type="response"), tolerance=.05, check.attributes=FALSE)

############################################
.test = "lasso/scad/mcp all seem to work" ##
############################################
n <- 100
p <- 10
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
X <- matrix(rnorm(n*p), n, p)
par(mfrow=c(3,1))

fit <- ncvsurv(X, y, penalty="lasso")
plot(fit, main="lasso")
fit <- ncvsurv(X, y, penalty="SCAD")
plot(fit, main="scad")
fit <- ncvsurv(X, y, penalty="MCP")
plot(fit, main="mcp")

################################
.test = "logLik() is correct" ##
################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
y <- Surv(rexp(n), rbinom(n, 1, 0.5))

fit.mle <- coxph(y~X)
fit <- ncvsurv(X, y, lambda.min=0)
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)

nfit <- ncvsurv(X, y, penalty="lasso")
l <- nfit$lambda
gfit <- glmnet(X, y, family="cox", lambda=l)
-2*logLik(nfit)[50]
glmnet:::coxnet.deviance(predict(gfit, X, s=l[50]), y)

#############################################
.test = "ncvsurv handles constant columns" ##
#############################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
X[, 3:5] <- 0
y <- Surv(rexp(n), rbinom(n, 1, 0.5))
fit <- ncvsurv(X, y)

#######################################
.test = "cv.ncvsurv() seems to work" ##
#######################################
n <- 50
p <- 100
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
par(mfrow=c(2,2))

ncvfit <- cv.ncvsurv(X, y, penalty="lasso", lambda.min=0.2)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, penalty="lasso", events.only=FALSE, lambda.min=0.2)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, lambda.min=0.4)
plot(ncvfit)
ncvfit <- cv.ncvsurv(X, y, events.only=FALSE, lambda.min=0.4)
plot(ncvfit)

require(glmnet)
par(mfrow=c(2,1))
ncvfit <- cv.ncvsurv(X, y, penalty="lasso")
plot(ncvfit)
gcvfit <- cv.glmnet(X, y, family="cox", grouped=TRUE, lambda=ncvfit$lambda)
plot(gcvfit)

#########################################
.test = "penalty.factor seems to work" ##
#########################################
n <- 50
p <- 4
X <- matrix(rnorm(n*p), ncol=p)
X <- prcomp(X)$x
y <- Surv(rexp(n, 1), rbinom(n, 1, 0.5))
penalty.factor=c(0,0,1,10)

par(mfrow=c(2,1))
fit <- ncvsurv(X, y)
plot(fit)
fit <- ncvsurv(X, y, penalty.factor=penalty.factor)
plot(fit)

######################################
.test = "cv.ncvsurv() options work" ##
######################################
n <- 50
p <- 20
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))

par(mfrow=c(2,2))
cvfit <- cv.ncvsurv(X, y)
plot(cvfit, type="all")
par(mfrow=c(2,2))
cvfit <- cv.ncvsurv(X, y, events.only=FALSE)
plot(cvfit, type="all")
print(summary(cvfit))
print(predict(cvfit, "coefficients"))
print(predict(cvfit, "vars"))
print(predict(cvfit, "nvars"))

#############################################################
.test = "cox survival predictions agree with Kaplan-Meier" ##
#############################################################
n <- 50
p <- 5
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, rep(0, p-2))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.75))

fit <- ncvsurv(X, y)
S <- predict(fit, X[1,], which=1, type='survival')
km <- survfit(y~1)
par(op)
plot(km, conf.int=FALSE, mark.time=FALSE, xlim=c(0,10), lwd=10, col="gray")
lines(fit$time, S(fit$time), type="s", col="slateblue", lwd=2)
median(km)
predict(fit, X[1,], which=1, type='median')

######################################################
.test = "cox survival predictions agree with coxph" ##
######################################################
n <- 50
p <- 5
# Standardize X, otherwise coxph makes strange adjustments for the mean
#X <- .Call('standardize', rbind(rep(0,p), matrix(rnorm((n-1)*p), ncol=p)))[[1]]
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, rep(0, p-2))
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
df <- data.frame(y, X)

fit <- ncvsurv(X, y, lambda.min=0)
sfit <- coxph(y~., data=df)
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (i in 1:9) {
  S <- predict(fit, X[i,], which=100, type='survival')
  km.cox <- survfit(sfit, newdata = df[i,], type="kalbfleisch-prentice")
  plot(km.cox, conf.int=FALSE, mark.time=FALSE, xlim=c(0,10), lwd=10, col="gray")
  lines(fit$time, S(fit$time), type="s", col="slateblue", lwd=2)
}

##########################################
.test = "cox survival predictions work" ##
##########################################
data(Lung, package='ncvreg')
X <- Lung$X
y <- Lung$y

fit <- ncvsurv(X, y)
M <- predict(fit, X, type='median')
M[1:10, 1:10]
S <- predict(fit, X, lambda=0.3, type='survival')
sapply(S, function(f) f(100))
par(op)
plot(S, xlim=c(0,200))
S <- predict(fit, X, lambda=0.05, type='survival')
plot(S, xlim=c(0,200))
S <- predict(fit, X, lambda=0.44, type='survival')
plot(S, xlim=c(0,200))

##################################
.test = "AUC calculation works" ##
##################################
data(Lung, package='ncvreg')
X <- Lung$X
y <- Lung$y

cvfit <- cv.ncvsurv(X, y, returnY=TRUE)
AUC <- auc(cvfit)
AUC[length(AUC)]
fit <- coxph(y~X[,1:8])
survConcordance(y~predict(fit))$concordance
