refresh(ncvreg)
X <- matrix(rnorm(500),ncol=10)
b <- rnorm(10)
y <- rnorm(X%*%b)
yy <- y > 0
tol <- .01
fit <- ncvreg(X,y)
fit <- ncvreg(X,yy,family="binomial")

coef <- lm(y~X)$coef
beta <- ncvreg(X,y,lambda=0,penalty="SCAD",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")
beta <- ncvreg(X,y,lambda=0,penalty="MCP",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")

coef <- glm(yy~X,family="binomial")$coef
beta <- ncvreg(X,yy,family="binomial",lambda=0,penalty="SCAD",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")
beta <- ncvreg(X,yy,family="binomial",lambda=0,penalty="MCP",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")

fit <- ncvreg(X,y)
fit <- ncvreg(X,yy,family="binomial")
