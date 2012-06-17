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
