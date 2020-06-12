# Linear
data(Prostate, package='ncvreg')
X <- Prostate$X
y <- Prostate$y
cvfit <- cv.ncvreg(X, y)
summary(cvfit)
plot(cvfit, type='rsq')
summary(lm(y~X))

# Logistic
data(Heart, package='ncvreg')
X <- Heart$X
y <- Heart$y
cvfit <- cv.ncvreg(X, y, family='binomial')
summary(cvfit)
plot(cvfit, type='rsq')
l1 = logLik(glm(y~X, family='binomial'))[1]
l0 = logLik(glm(y~1, family='binomial'))[1]
1 - exp(-2/length(y) * (l1 - l0))

# Cox
data(Lung, package='ncvreg')
X <- Lung$X
y <- Lung$y
cvfit <- cv.ncvsurv(X, y)
summary(cvfit)
plot(cvfit, type='rsq')
library(survival)
summary(coxph(y~X))
