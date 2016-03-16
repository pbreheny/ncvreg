#source("~/dev/.ncvreg.setup.R")

## Gen data
n <- 40
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
y <- rnorm(n, mean=X%*%b, sd=2)

####################################################
.test = "perm.ncvreg works for linear regression" ##
####################################################
pmfit <- perm.ncvreg(X, y)

par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)

######################################################
.test = "perm.ncvreg works for logistic regression" ##
######################################################
pmfit <- perm.ncvreg(X, y > 0, family='binomial')

par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)
plot(pmfit$fit, log=TRUE)

#####################################################
.test = "perm.ncvreg works for Poisson regression" ##
#####################################################
pmfit <- perm.ncvreg(X, rank(y), family='poisson')

par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)
plot(pmfit$fit, log=TRUE)

##########################
.test = "permres works" ##
##########################
fit <- ncvreg(X, y, returnX=TRUE)
permres.ncvreg(fit, lam=0.2)

#############################################################
.test = "permute='residuals' option for perm.ncvreg works" ##
#############################################################
pmfit <- perm.ncvreg(X, y, permute="residuals", N=25)
par(mfrow=c(2,2))
plot(pmfit)
plot(pmfit, type="EF")
plot(pmfit$fit)
plot(pmfit$fit)
