
library("coin")
set.seed(290875)

cort <- function(n = 100, ties = FALSE) {
    if (ties) {
        mydf <- data.frame(x = runif(n), y = round(runif(n), 1))
        dist <- asymptotic()
    } else {
        mydf <- data.frame(x = runif(n), y = runif(n))
        dist <- approximate(B = 9999)
    }
    p1 <- cor.test(~ x + y, data = mydf, method = "spearman")$p.value
    p2 <- pvalue(spearman_test(x ~ y, data = mydf, dist = dist))
    p1 - p2
}

p <- numeric(100)
for (i in 1:length(p))
    p[i] <- cort()
summary(abs(p))

p <- numeric(100)
for (i in 1:length(p))
    p[i] <- cort(ties = TRUE)
summary(abs(p))

