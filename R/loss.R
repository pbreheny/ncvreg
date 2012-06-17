loss.ncvreg <- function(y,yhat,family)
  {
    n <- length(y)
    if (family=="gaussian") return(apply(y-yhat,2,crossprod)/(2*n))
    if (family=="binomial") return((apply(log(yhat[y==1,,drop=FALSE]),2,sum)+apply(log(1-yhat[y==0,,drop=FALSE]),2,sum))/(-n))
  }
