source("~/dev/.ncvreg.setup.R")

test_that("ncvreg works for linear regression", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), n, p)
  b <- rnorm(p)
  y <- rnorm(n, X%*%b)
  beta <- lm(y~X)$coef
  scad <- coef(ncvreg(X,y,lambda=0,penalty="SCAD",eps=.0001))
  mcp <- coef(ncvreg(X,y,lambda=0,penalty="MCP",eps=.0001))
  expect_that(scad,equals(beta,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(beta,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg works for logistic regression", {
  X <- matrix(rnorm(500),ncol=10)
  b <- rnorm(10)
  y <- rnorm(X%*%b) > 0
  beta <- glm(y~X,family="binomial")$coef
  scad <- coef(ncvreg(X,y,lambda=0,family="binomial",penalty="SCAD",eps=.0001))
  mcp <- coef(ncvreg(X,y,lambda=0, family="binomial",penalty="MCP", eps=.0001))
  expect_that(scad,equals(beta,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(beta,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg works for Poisson regression", {
  X <- matrix(rnorm(500),ncol=10)
  b <- rnorm(10)
  y <- rnorm(X%*%b) > 0
  beta <- glm(y~X,family="poisson")$coef
  scad <- coef(ncvreg(X,y,lambda=0,family="poisson",penalty="SCAD",eps=.0001))
  mcp <- coef(ncvreg(X,y,lambda=0, family="poisson",penalty="MCP", eps=.0001))
  expect_that(scad,equals(beta,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(beta,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg reproduces lasso", {
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
  expect_that(nlasso, equals(glasso, tolerance=.01, check.attributes=FALSE))
  
  yy <- runif(n) > .5
  nlasso <- coef(fit <- ncvreg(X, yy, family="binomial", penalty="lasso"))
  plot(fit, log=TRUE)
  glasso <- as.matrix(coef(fit <- glmnet(X, yy, family="binomial", lambda=fit$lambda)))
  plot(fit, "lambda")  
  expect_that(nlasso, equals(glasso, tolerance=.01, check.attributes=FALSE))

  yy <- rpois(n, 1)
  nlasso <- coef(fit <- ncvreg(X, yy, family="poisson", penalty="lasso"))
  plot(fit, log=TRUE)
  glasso <- as.matrix(coef(fit <- glmnet(X, yy, family="poisson", lambda=fit$lambda)))
  plot(fit, "lambda")  
  expect_that(nlasso, equals(glasso, tolerance=.01, check.attributes=FALSE))
})

test_that("logLik() is correct", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  
  y <- rnorm(n)
  fit.mle <- lm(y~X)
  fit <- ncvreg(X, y, lambda.min=0)
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, tol=.001))
  
  yy <- runif(n) > .5
  fit.mle <- glm(yy~X, family="binomial")
  fit <- ncvreg(X, yy, lambda.min=0, family="binomial")
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, , tol=.001))
  
  yy <- rpois(n, 1)
  fit.mle <- glm(yy~X, family="poisson")
  fit <- ncvreg(X, yy, lambda.min=0, family="poisson")
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(logLik(fit))[100], equals(AIC(fit.mle), check.attributes=FALSE, tol=.001))
})

test_that("ncvreg handles constant columns", {
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
})

test_that("cv.ncvreg() seems to work", {
  n <- 50
  p <- 100
  X <- matrix(rnorm(n*p), ncol=p)
  b <- c(2, -2, 1, -1, rep(0, p-4))
  y <- rnorm(n, mean=X%*%b, sd=2)
  yb <- y > .5
  yp <- rpois(n, exp(X%*%b/3))
  par(mfrow=c(3,2))

  require(glmnet)
  cvfit <- cv.glmnet(X, y)
  plot(cvfit)
  cvfit <- cv.ncvreg(X, y, penalty="lasso")
  plot(cvfit)
  cvfit <- cv.glmnet(X, yb, family="binomial", lambda.min=0)
  plot(cvfit)
  cvfit <- cv.ncvreg(X, yb, family="binomial", penalty="lasso", lambda.min=0)
  plot(cvfit)
  cvfit <- cv.glmnet(X, yp, family="poisson", lambda.min=0)
  plot(cvfit)
  cvfit <- cv.ncvreg(X, yp, family="poisson", penalty="lasso", lambda.min=0)
  plot(cvfit)
})  

test_that("penalty.factor seems to work", {
  n <- 50
  p <- 4
  X <- matrix(rnorm(n*p), ncol=p)
  y <- rnorm(n)
  yy <- y > .5
  penalty.factor=c(0,0,1,10)
  
  par(mfrow=c(2,2))
  fit <- ncvreg(X, y)
  plot(fit)
  fit <- ncvreg(X, y, penalty.factor=penalty.factor)
  plot(fit)
  fit <- ncvreg(X, yy, family="binomial")
  plot(fit)
  fit <- ncvreg(X, yy, family="binomial", penalty.factor=penalty.factor)
  plot(fit)
})

test_that("cv.ncvreg() options work for gaussian", {
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
})  

test_that("cv.ncvreg() options work for binomial", {
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
})  

test_that("cv.ncvreg() options work for poisson", {
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
})  

test_that("ncvreg_fit works", {
  n <- 200
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  b <- c(-3, 3, rep(0, 8))
  y <- rnorm(n, mean=X%*%b, sd=1)
  yy <- y > 0
  
  b1 <- ncvreg_fit(X, y, penalty="lasso", lam=c(1, 0.1, 0.01))
  b2 <- glmnet(X, y, lambda=c(1, 0.1, 0.01), standardize=FALSE, intercept=FALSE)
  expect_that(b1$beta, equals(as.matrix(coef(b2))[-1,], check.attributes=FALSE, tol=.001))
  
  b1 <- ncvreg_fit(X, y, "binomial", "lasso", lam=c(1, 0.1, 0.01))
  b2 <- glmnet(X, y, "binomial", lambda=c(1, 0.1, 0.01), standardize=FALSE, intercept=TRUE)
  expect_that(b1$beta, equals(as.matrix(coef(b2)), check.attributes=FALSE, tol=.001))
  
  b1 <- ncvreg_fit(X, y, "poisson", "lasso", lam=c(1, 0.1, 0.01))
  b2 <- glmnet(X, y, "poisson", lambda=c(1, 0.1, 0.01), standardize=FALSE, intercept=TRUE)
  expect_that(b1$beta, equals(as.matrix(coef(b2)), check.attributes=FALSE, tol=.001))

  b1 <- ncvreg_fit(X, y, penalty="MCP", lam=0.01)
  b2 <- ncvreg_fit(X, y, "binomial", penalty="MCP", lam=c(1, 0.1, 0.01))
  b3 <- ncvreg_fit(X, y, penalty="SCAD", lam=c(1, 0.1, 0.01))
  b4 <- ncvreg_fit(X, y, "binomial", penalty="SCAD", lam=c(1, 0.1, 0.01))
})

