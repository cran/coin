
library("coin")

set.seed(290875)

### cmh 2 x 2 x k

cmh22k <- function(k = 5, n = 100, alternative = "two.sided") {

    mydf <- data.frame(x = gl(2, n / 2)[sample(1:n)],
                       y = gl(2, n / 2)[sample(1:n)],
                       z = gl(k, n / k))
    p1 <- mantelhaen.test(mydf$x, mydf$y, mydf$z, correct = FALSE,
                          alternative = alternative, exact = TRUE)$p.value
    p2 <- pvalue(independence_test(y ~ x | z, data = mydf, dist = approximate(B = 29999), 
                 alternative = alternative))
    p1 - p2
}

k <- rep(c(2, 4, 5, 10), 2)
p <- numeric(100)

### two-sided
for (i in 1:length(p)) {
    p[i] <- cmh22k(k = sample(k)[1], alt = "t", n = 100)
    print(p[i])
}
summary(abs(p))

### less
for (i in 1:length(p)) {
    p[i] <- cmh22k(k = sample(k)[1], alt = "l", n = 100)
    print(p[i])
}
summary(abs(p))

### greater
for (i in 1:length(p)) {
    p[i] <- cmh22k(k = sample(k)[1], alt = "g", n = 100)
    print(p[i])
}
summary(abs(p))


### cmh i x j x k
cmhijk <- function(k = c(2, 2, 5), n = 100, alternative = "two.sided") {

    mydf <- data.frame(x = gl(k[1], n / k[1])[sample(1:n)],
                       y = gl(k[2], n / k[2])[sample(1:n)],
                       z = gl(k[3], n / k[3]))
    p1 <- mantelhaen.test(mydf$x, mydf$y, mydf$z, correct = FALSE,
                          alternative = alternative, exact = FALSE)$p.value
    p2 <- pvalue(cmh_test(y ~ x | z, data = mydf, alternative = alternative))
    p1 - p2
}

### two-sided
k <- rep(c(2, 4, 5, 10), 2)
p <- numeric(100)
for (i in 1:length(p)) {
    p[i] <- cmhijk(sample(k)[1:3], alt = "t", n = 100)
    print(p[i])
}
summary(abs(p))


