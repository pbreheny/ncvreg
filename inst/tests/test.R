source("~/dev/.ncvreg.setup.R")

test_that("ncvreg works for linear regression", {
  X <- matrix(rnorm(500),ncol=10)
  b <- rnorm(10)
  y <- rnorm(X%*%b)
  coef <- lm(y~X)$coef
  scad <- coef(ncvreg(X,y,lambda=0,penalty="SCAD",eps=.0001))
  mcp <- coef(ncvreg(X,y,lambda=0,penalty="MCP",eps=.0001))
  expect_that(scad,equals(coef,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(coef,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg works for logistic regression", {
  X <- matrix(rnorm(500),ncol=10)
  b <- rnorm(10)
  y <- rnorm(X%*%b) > 0
  coef <- glm(y~X,family="binomial")$coef
  scad <- coef(ncvreg(X,y,lambda=0,family="binomial",penalty="SCAD",eps=.0001))
  mcp <- coef(ncvreg(X,y,lambda=0, family="binomial",penalty="MCP", eps=.0001))
  expect_that(scad,equals(coef,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(coef,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg reproduces lasso", {
  require(glmnet)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  y <- rnorm(n)
  yy <- runif(n) > .5
  nlasso <- coef(fit <- ncvreg(X, y, gamma=1e8))
  par(mfrow=c(2,2)); plot(fit, log=TRUE)
  glasso <- as.matrix(coef(fit <- glmnet(X, y, lambda=fit$lambda)))
  plot(fit, "lambda")
  expect_that(nlasso, equals(glasso, tolerance=.01, check.attributes=FALSE))
  nlasso <- coef(fit <- ncvreg(X, yy, family="binomial", gamma=1e8))
  plot(fit, log=TRUE)
  glasso <- as.matrix(coef(fit <- glmnet(X, yy, family="binomial", lambda=fit$lambda)))
  plot(fit, "lambda")  
  expect_that(nlasso, equals(glasso, tolerance=.01, check.attributes=FALSE))
})

test_that("logLik() is correct", {
  n <- 50
  p <- length(10)
  X <- matrix(rnorm(n*p), ncol=p)
  y <- rnorm(n)
  yy <- runif(n) > .5
  fit.mle <- lm(y~X)
  fit <- ncvreg(X, y, lambda.min=0)
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, tol=.001))
  fit.mle <- glm(yy~X, family="binomial")
  fit <- ncvreg(X, yy, lambda.min=0, family="binomial")
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, , tol=.001))
})

test_that("ncvreg handles constant columns", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  y <- rnorm(n)
  yy <- runif(n) > .5
  X[, 3:5] <- 0
  fit <- ncvreg(X, y)
  fit <- ncvreg(X, yy, family="binomial")
})

test_that("cv.ncvreg() seems to work", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  b <- rnorm(p, sd=2)
  b[abs(b) < 1] <- 0
  y <- rnorm(n, mean=X%*%b)
  yy <- y > .5
  
  par(mfrow=c(2,2))
  require(glmnet)
  cvfit <- cv.glmnet(X, y)
  plot(cvfit)
  cvfit <- cv.ncvreg(X, y, gamma=1e8)
  plot(cvfit)
  cvfit <- cv.glmnet(X, yy, family="binomial", lambda.min=0)
  plot(cvfit)
  cvfit <- cv.ncvreg(X, yy, family="binomial", gamma=1e8, lambda.min=0)
  plot(cvfit)
})  

## plot
## cv.ncvreg

