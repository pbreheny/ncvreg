require(survival)
set.seed(1)
equal <- function(x, y) {all.equal(x, y, tol=0.001, check.attributes=FALSE)}

# Works
y <- Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*10), 50, 10)
fit <- ncvsurv(X, y, lambda.min=0)

# Equals MLE when lam=0
fit.mle <- coxph(y~X)
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# logLik
stopifnot(equal(logLik(fit)[100], logLik(fit.mle)[1]))
stopifnot(equal(AIC(fit)[100], AIC(fit.mle)))

# Other penalties
fit <- ncvsurv(X, y, lambda.min=0, penalty="SCAD")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))
fit <- ncvsurv(X, y, lambda.min=0, penalty="lasso")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# Predict
p <- predict(fit, X, 'link', lambda=0.1)
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')

# Penalty factor
fit <- ncvsurv(X, y, penalty.factor=c(0:9))

if (FALSE) {
  #############################################################
  .test = "cox survival predictions agree with Kaplan-Meier" ##
  #############################################################
  n <- 50
  p <- 5
  X <- matrix(rnorm(n*p), ncol=p)
  b <- c(2, -2, rep(0, p-2))
  y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.75))

  fit <- ncvsurv(X, y)
  S <- predict(fit, X[1,], which=1, type='survival')
  km <- survfit(y~1)
  par(op)
  plot(km, conf.int=FALSE, mark.time=FALSE, xlim=c(0,10), lwd=10, col="gray")
  lines(fit$time, S(fit$time), type="s", col="slateblue", lwd=2)
  median(km)
  predict(fit, X[1,], which=1, type='median')

  ######################################################
  .test = "cox survival predictions agree with coxph" ##
  ######################################################
  n <- 50
  p <- 5
  # Standardize X, otherwise coxph makes strange adjustments for the mean
  #X <- .Call('standardize', rbind(rep(0,p), matrix(rnorm((n-1)*p), ncol=p)))[[1]]
  X <- matrix(rnorm(n*p), ncol=p)
  b <- c(2, -2, rep(0, p-2))
  y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.5))
  df <- data.frame(y, X)

  fit <- ncvsurv(X, y, lambda.min=0)
  sfit <- coxph(y~., data=df)
  par(mfrow=c(3,3), mar=c(1,1,1,1))
  for (i in 1:9) {
    S <- predict(fit, X[i,], which=100, type='survival')
    km.cox <- survfit(sfit, newdata = df[i,], type="kalbfleisch-prentice")
    plot(km.cox, conf.int=FALSE, mark.time=FALSE, xlim=c(0,10), lwd=10, col="gray")
    lines(fit$time, S(fit$time), type="s", col="slateblue", lwd=2)
  }

  ##########################################
  .test = "cox survival predictions work" ##
  ##########################################
  data(Lung, package='ncvreg')
  X <- Lung$X
  y <- Lung$y

  fit <- ncvsurv(X, y)
  M <- predict(fit, X, type='median')
  M[1:10, 1:10]
  S <- predict(fit, X, lambda=0.3, type='survival')
  sapply(S, function(f) f(100))
  par(op)
  plot(S, xlim=c(0,200))
  S <- predict(fit, X, lambda=0.05, type='survival')
  plot(S, xlim=c(0,200))
  S <- predict(fit, X, lambda=0.44, type='survival')
  plot(S, xlim=c(0,200))
}
