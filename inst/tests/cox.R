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
  
  nlasso <- coef(fit <- ncvreg(X, y, penalty="lasso", family="cox"))
  ##nlasso <- coef(fit <- ncvreg(X, y, penalty="lasso", family="cox", lambda=fit$lambda))
  plot(fit, log=TRUE)
  glasso <- as.matrix(coef(fit <- glmnet(X, y, family="cox", lambda=fit$lambda)))
  plot(fit, "lambda")
  expect_that(nlasso, equals(glasso, tolerance=.01, check.attributes=FALSE))  
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
