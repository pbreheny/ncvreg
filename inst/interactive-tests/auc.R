#source("~/dev/.ncvreg.setup.R")
require(ROCR)

##############################################
.test = "cv.ncvreg() agrees with AUC" ##
##############################################
n <- 100
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-1, 1, rep(0, 8))

y <- rnorm(n, mean=X%*%b) > 0
cvfit <- cv.ncvreg(X, y, family='binomial', returnY=TRUE)
L <- ncol(cvfit$Y)
auc <- numeric(L)
for (j in 1:L) {
  pred <- prediction(cvfit$Y[,j], y)
  auc[j] <- performance(pred, 'auc')@y.values[[1]]
}
ll <- log(cvfit$lambda)
plot(ll, auc, xlim=rev(range(ll)), type='l', lwd=3)
calcAUC <- function(y, P) {
  if (!('dim' %in% names(attributes(P)))) P <- matrix(P, ncol=1)
  n <- sum(y)
  m <- length(y) - n
  J <- ncol(P)
  auc <- numeric(J)
  for (j in 1:J) {
    W <- sum(rank(P[,j])[y==1])
    auc[j] <- (W - n*(n+1)/2)/(m*n)
  }
  auc
}
calcAUC(y, cvfit$Y)
auc

###################################
.test = "cv.ncvsurv() returns Y" ##
###################################
n <- 50
p <- 20
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, 1, -1, rep(0, p-4))
#b <- rep(0, p)
y <- cbind(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))

cvfit <- cv.ncvsurv(X, y, returnY=TRUE, nfolds=n)
plot(cvfit)
fit <- ncvsurv(X[-1,], y[-1,], lambda=cvfit$lambda)
head(predict(fit, X[1,,drop=FALSE]))
head(cvfit$Y[1,])

ll <- log(cvfit$lambda)
cve <- apply(cvfit$Y, 2, function(x) loss.ncvsurv(y, x))
plot(ll, cve, xlim=rev(range(ll)), type='l', lwd=3, las=1)
