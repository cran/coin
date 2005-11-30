
library("coin")
set.seed(290875)

chisim <- function(k = c(5, 4), n = 100) {
    x <- gl(k[1], n / k[1])
    y <-  gl(k[2], n/k[2])[sample(1:n)]
    tab <- table(x, y)              
    p1 <- chisq.test(tab, correct = FALSE, simulate = TRUE, B = 9999)$p.value
    p2 <- pvalue(chisq_test(tab, dist = approximate(B = 9999)))
    if (abs(p1 - p2) > 0.1) save(tab, x, y, file = "tab.rda")
    p1 - p2 
}

p <- numeric(10)
k <- rep(c(2, 4, 5, 10), 2)
for (i in 1:length(p))
    p[i] <- chisim(k = sample(k)[1:2])

summary(abs(p))

chia <- function(k = c(5, 4), n = 100) {
    tab <- table(gl(k[1], n / k[1]), gl(k[2], n/k[2])[sample(1:n)])              
    p1 <- chisq.test(tab, correct = FALSE)$p.value
    p2 <- pvalue(chisq_test(tab))
    p1 - p2 
}

k <- rep(c(2, 4, 5, 10), 2)
for (i in 1:length(p))
    p[i] <- chia(k = sample(k)[1:2])

summary(abs(p))


