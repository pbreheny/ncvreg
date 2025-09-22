# X data frame, y factor
data(Heart, package = "ncvreg")
x <- Heart$X
y <- factor(Heart$y, labels = c("No", "Yes"))
fit <- ncvreg(x, y, family = "binomial")
s <- std(x)

# integer X, y
x <- as.matrix(round(x))
storage.mode(x) <- "integer"
fit <- ncvreg(x, Heart$y, family = "binomial")

# logical y
y <- Heart$y == 1
fit <- ncvreg(x, y, family = "binomial")
