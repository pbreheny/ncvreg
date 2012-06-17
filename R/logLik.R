logLik.ncvreg <- function(object,REML=TRUE)
  {
    n <- as.numeric(object$n)
    df <- as.numeric(apply(object$beta!=0,2,sum))
    if (object$family=="gaussian")
      {
        if (REML) rdf <- n-df
        else rdf <- n
        RSS <- object$loss
        l <- -n/2 * (log(2*pi) + log(RSS) - log(rdf)) - rdf/2
        df <- df + 1
      }
    if (object$family=="binomial") l <- -1*object$loss

    val <- l
    attr(val,"df") <- df
    attr(val,"nobs") <- n
    class(val) <- "logLik"
    return(val)    
  }
