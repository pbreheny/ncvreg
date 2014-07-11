source("~/dev/.ncvreg.setup.R")
require(survival)

test_that("ncvreg works for simple cox regression, no censoring", {
  n <- 10
  X <- matrix(runif(n), n, 1)
  y <- Surv(sort(rexp(n, X), decreasing=TRUE), rep(1,n))
  X <- X - mean(X)
  beta <- coef(coxph(y ~ X))
  eta <- X%*%beta
  w <- exp(eta)/cumsum(exp(eta))
  sum(X * (1-w))
  scad <- coef(ncvreg(X, y, family="cox", lambda=0, penalty="SCAD", eps=.0001))
  mcp <- coef(ncvreg(X, y, family="cox", lambda=0, penalty="MCP", eps=.0001))
  expect_that(scad,equals(beta,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(beta,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg works for simple cox regression", {
  n <- 30
  y <- Surv(rexp(n), rbinom(n, 1, 0.5))
  y[which.max(y[,1]), 2] <- 0
  y[which.min(y[,1]), 2] <- 0
  X <- matrix(rnorm(n), n, 1)
  beta <- coef(coxph(y ~ X))
  scad <- coef(ncvreg(X, y, family="cox", lambda=0, penalty="SCAD", eps=.0001))
  mcp <- coef(ncvreg(X, y, family="cox", lambda=0,penalty="MCP",eps=.0001))
  expect_that(scad,equals(beta,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(beta,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg works for cox regression", {
  n <- 100
  p <- 10
  y <- Surv(rexp(n), rbinom(n, 1, 0.5))
  X <- matrix(rnorm(n*p), n, p)
  beta <- coef(coxph(y ~ X))
  scad <- coef(ncvreg(X, y, family="cox", lambda=0, penalty="SCAD", eps=.0001))
  mcp <- coef(ncvreg(X, y, family="cox", lambda=0,penalty="MCP",eps=.0001))
  expect_that(scad,equals(beta,tolerance=.01,check.attributes=FALSE))
  expect_that(mcp,equals(beta,tolerance=.01,check.attributes=FALSE))
})

test_that("ncvreg agrees with coxnet", {
  require(glmnet)
  n <- 100
  p <- 10
  y <- Surv(rexp(n), rbinom(n, 1, 0.5))
  X <- matrix(rnorm(n*p), n, p)
  par(mfrow=c(2,1))
  
  nlasso <- coef(nfit <- ncvreg(X, y, penalty="lasso", family="cox"))
  ##nlasso <- coef(fit <- ncvreg(X, y, penalty="lasso", family="cox", lambda=fit$lambda))
  plot(nfit, log=TRUE)
  glasso <- as.matrix(coef(gfit <- glmnet(X, y, family="cox", lambda=nfit$lambda)))
  plot(gfit, "lambda")
  expect_that(nlasso, equals(glasso, tolerance=.01, check.attributes=FALSE))
  
  expect_that(predict(nfit, X, "link"), equals(predict(gfit, X, type="link"), tolerance=.01, check.attributes=FALSE))
  expect_that(predict(nfit, X, "response"), equals(predict(gfit, X, type="response"), tolerance=.01, check.attributes=FALSE))
})

test_that("lasso/scad/mcp", {
  n <- 100
  p <- 10
  y <- Surv(rexp(n), rbinom(n, 1, 0.5))
  X <- matrix(rnorm(n*p), n, p)
  par(mfrow=c(3,1))
  
  fit <- ncvreg(X, y, penalty="lasso", family="cox")
  plot(fit, main="lasso")
  fit <- ncvreg(X, y, penalty="SCAD", family="cox")
  plot(fit, main="scad")
  fit <- ncvreg(X, y, penalty="MCP", family="cox")
  plot(fit, main="mcp")
})

test_that("logLik() is correct", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  y <- Surv(rexp(n), rbinom(n, 1, 0.5))
  
  fit.mle <- coxph(y~X)
  fit <- ncvreg(X, y, lambda.min=0, family="cox")
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, tol=.001))
})

test_that("ncvreg handles constant columns", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  X[, 3:5] <- 0
  y <- Surv(rexp(n), rbinom(n, 1, 0.5))
  fit <- ncvreg(X, y, family="cox")
})

test_that("cv.ncvreg() seems to work", {
  n <- 50
  p <- 100
  X <- matrix(rnorm(n*p), ncol=p)
  b <- c(2, -2, 1, -1, rep(0, p-4))
  y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
  par(mfrow=c(2,1))
  
  require(glmnet)
  cvfit <- cv.glmnet(X, y, family="cox")
  plot(cvfit)
  cvfit <- cv.ncvreg(X, y, family="cox", penalty="lasso")
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

test_that("ncvreg_fit works", {
  n <- 200
  p <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  b <- c(-3, 3, rep(0, 8))
  y <- rnorm(n, mean=X%*%b, sd=1)
  
  b1 <- ncvreg_fit(X, y, penalty="lasso", lam=c(1, 0.1, 0.01))
  b2 <- glmnet(X, y, lambda=c(1, 0.1, 0.01), standardize=FALSE, intercept=FALSE)
  expect_that(b1$beta, equals(as.matrix(coef(b2))[-1,], check.attributes=FALSE, tol=.001))
  
  yy <- y > 0
  b1 <- ncvreg_fit(X, yy, "binomial", "lasso", lam=c(1, 0.1, 0.01))
  b2 <- glmnet(X, yy, "binomial", lambda=c(1, 0.1, 0.01), standardize=FALSE, intercept=TRUE)
  expect_that(b1$beta, equals(as.matrix(coef(b2)), check.attributes=FALSE, tol=.001))
  
  b1 <- ncvreg_fit(X, yy, "poisson", "lasso", lam=c(1, 0.1, 0.01))
  b2 <- glmnet(X, yy, "poisson", lambda=c(1, 0.1, 0.01), standardize=FALSE, intercept=TRUE)
  expect_that(b1$beta, equals(as.matrix(coef(b2)), check.attributes=FALSE, tol=.001))
  
  b1 <- ncvreg_fit(X, y, penalty="MCP", lam=0.01)
  b2 <- ncvreg_fit(X, yy, "binomial", penalty="MCP", lam=c(1, 0.1, 0.01))
  b3 <- ncvreg_fit(X, y, penalty="SCAD", lam=c(1, 0.1, 0.01))
  b4 <- ncvreg_fit(X, yy, "binomial", penalty="SCAD", lam=c(1, 0.1, 0.01))
})
