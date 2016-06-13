data(heart)
X <- heart[,1:9]
y <- factor(heart$chd, labels=c("No", "Yes"))
fit <- ncvreg(X, y, family="binomial")
std(X)
