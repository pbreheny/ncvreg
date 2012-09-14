setupLambda <- function(X, y, family, alpha, lambda.min, nlambda)
  {
    n <- nrow(X)
    p <- ncol(X)

    ## Determine lambda.max
    if (family=="gaussian")
      {
        r <- y - mean(y)
        l1.max <- max(abs(crossprod(X,r)/n))    
      }
    if (family=="binomial")
      {
        fit <- glm(y~1,family="binomial")
        w <- fit$weights
        r <- residuals(fit, "working")
        l1.max <- max(abs(crossprod(X,w*r)/n))
      }
    lambda.max <- l1.max/alpha
    
    if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
    else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
    return(lambda)
  }
