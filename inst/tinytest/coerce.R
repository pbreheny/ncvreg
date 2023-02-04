# X data frame, y factor
data(Heart, package='ncvreg')
X <- Heart$X
y <- factor(Heart$y, labels=c("No", "Yes"))
fit <- ncvreg(X, y, family="binomial")
s <- std(X)

# integer X, y
X <- as.matrix(round(X))
storage.mode(X) <- "integer"
fit <- ncvreg(X, Heart$y, family="binomial")

# logical y
y <- Heart$y == 1
fit <- ncvreg(X, y, family="binomial")
