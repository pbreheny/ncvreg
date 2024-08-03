# Works

## lasso
X <- matrix(rnorm(500), 50, 10)
y <- X[,1] + rnorm(50)
boot <- boot_ncvreg(X, y, alpha = 0.5)

## MCP
boot <- boot_ncvreg(X, y, penalty = "MCP")

## SCAD
boot <- boot_ncvreg(X, y, penalty = "SCAD")

## Elastic Net
boot <- boot_ncvreg(X, y, alpha = 0.5)