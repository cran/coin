
library("libcoin")
library("coin")
set.seed(29)

n <- 100
p <- 4
q <- 2
X <- matrix(runif(p * n), nc = p)
Y <- matrix(runif(q * n), nc = q)
w <- as.integer(floor(runif(n, max = 4)))
s <- sample(1:n, floor(n/2), replace = TRUE)
b <- sample(gl(2, 2, length = n))

cmp <- function(t1, t2) {
    if (is.null(t1$Covariance)) {
        var1 <- t1$Variance
        var2 <- diag(covariance(t2))
    } else {
        var1 <- t1$Covariance
        var2 <- covariance(t2)
        var2 <- var2[!upper.tri(var2)]
    }
    all.equal(
        list(t1$LinearStatistic, t1$Expectation, var1),
        list(t2@statistic@linearstatistic, t2@statistic@expectation, var2),
        check.attributes = FALSE
    )
}

cmp2 <- function(t1, t2) {
    nm <- c("LinearStatistic", "Expectation",
            if(t1$varonly == 1) "Variance" else "Covariance")
    all.equal(t1[nm], t2[nm])
}

t1 <-LinStatExpCov(X, Y)
t1v <-LinStatExpCov(X, Y, varonly = TRUE)
t2 <- independence_test(Y ~ X)
cmp(t1, t2)
cmp(t1v, t2)

t1 <-LinStatExpCov(X, Y, weights = w)
t1v <-LinStatExpCov(X, Y, weights = w, varonly = TRUE)
t2 <- independence_test(Y ~ X, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, subset = s)
t1v <- LinStatExpCov(X, Y, subset = s, varonly = TRUE)
t2 <- independence_test(Y ~ X, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, weights = w, subset = s)
t1v <- LinStatExpCov(X, Y, weights = w, subset = s, varonly = TRUE)
t2 <- independence_test(Y ~ X, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, block = b)
t1v <- LinStatExpCov(X, Y, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X  | b)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, weights = w, block = b)
t1v <- LinStatExpCov(X, Y, weights = w, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, weights = w, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, weights = w, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

n <- 100
n1 <- 5
n2 <- 4
p <- 4
q <- 2
X <- rbind(0, matrix(runif(p * n1), nc = p))
Y <- rbind(0, matrix(runif(q * n2), nc = q))
ix <- sample(1:n1, n, replace = TRUE)
iy <- sample(1:n2, n, replace = TRUE)
w <- as.integer(floor(runif(n, max = 4)))
s <- sample(1:n, floor(n/2), replace = TRUE)
b <- sample(gl(2, 2, length = n))

t1 <- LinStatExpCov(X, Y, ix, iy)
t1v <- LinStatExpCov(X, Y, ix, iy, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,])
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, ix, iy, weights = w)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,], weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, ix, iy, subset = s)
t1v <- LinStatExpCov(X, Y, ix, iy, subset = s, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,], subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,], weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, ix, iy, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,]  | b)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, ix, iy, weights = w, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,] | b, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, ix, iy, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,] | b, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,]| b, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

#### X factor
n <- 10000
p <- 40
q <- 20
X <- diag(p)[fx <- unclass(factor(sample(1:p, n, replace = TRUE))),]
Y <- matrix(runif(q * n), nc = q)
w <- as.integer(floor(runif(n, max = 4)))
s <- sample(1:n, floor(n/2), replace = TRUE)
b <- sample(gl(2, 2, length = n))

t1 <-LinStatExpCov(X, Y)
t1v <-LinStatExpCov(X, Y, varonly = TRUE)
t2 <- independence_test(Y ~ X)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y)
t1vf <-LinStatExpCov(fx, Y, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <-LinStatExpCov(X, Y, weights = w)
t1v <-LinStatExpCov(X, Y, weights = w, varonly = TRUE)
t2 <- independence_test(Y ~ X, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y, weights = w)
t1vf <-LinStatExpCov(fx, Y, weights = w, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, subset = s)
t1v <- LinStatExpCov(X, Y, subset = s, varonly = TRUE)
t2 <- independence_test(Y ~ X, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y, subset = s)
t1vf <-LinStatExpCov(fx, Y, subset = s, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, weights = w, subset = s)
t1v <- LinStatExpCov(X, Y, weights = w, subset = s, varonly = TRUE)
t2 <- independence_test(Y ~ X, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y, weights = w, subset = s)
t1vf <-LinStatExpCov(fx, Y, weights = w, subset = s, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, block = b)
t1v <- LinStatExpCov(X, Y, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X  | b)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y, block = b)
t1vf <-LinStatExpCov(fx, Y, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)


t1 <- LinStatExpCov(X, Y, weights = w, block = b)
t1v <- LinStatExpCov(X, Y, weights = w, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y, weights = w, block = b)
t1vf <-LinStatExpCov(fx, Y, weights = w, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y, subset = s, block = b)
t1vf <-LinStatExpCov(fx, Y, subset = s, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, weights = w, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, weights = w, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(fx, Y, weights = w, subset = s, block = b)
t1vf <-LinStatExpCov(fx, Y, weights = w, subset = s, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

n1 <- 5
n2 <- 7
X <- rbind(0, diag(n1))
Y <- rbind(0, matrix(runif(q * n2), nc = q))
ix <- sample(1:n1, n, replace = TRUE)
iy <- sample(1:n2, n, replace = TRUE)
w <- as.integer(floor(runif(n, max = 4)))
s <- sample(1:n, floor(n/2), replace = TRUE)
b <- sample(gl(2, 2, length = n))

YY <- Y[iy + 1,]
XX <- X[ix + 1,]

t1 <- LinStatExpCov(X, Y, ix, iy)
t1v <- LinStatExpCov(X, Y, ix, iy, varonly = TRUE)
t2 <- independence_test(YY ~ XX)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, ix, iy, weights = w)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, varonly = TRUE)
t2 <- independence_test(YY ~ XX, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, ix, iy, subset = s)
t1v <- LinStatExpCov(X, Y, ix, iy, subset = s, varonly = TRUE)
t2 <- independence_test(YY ~ XX, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy, subset = s)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, subset = s, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)


t1 <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s, varonly = TRUE)
t2 <- independence_test(YY ~ XX, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w, subset = s)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w, subset = s, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, ix, iy, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, block = b, varonly = TRUE)
t2 <- independence_test(YY ~ XX  | b)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy, block = b)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, ix, iy, weights = w, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, block = b, varonly = TRUE)
t2 <- independence_test(YY ~ XX | b, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w, block = b)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, ix, iy, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(YY ~ XX | b, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy, subset = s, block = b)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, subset = s, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

t1 <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, ix, iy, weights = w, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(YY ~ XX| b, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
t1f <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w, subset = s, block = b)
t1vf <-LinStatExpCov(numeric(0), Y, ix, iy, weights = w, subset = s, block = b, varonly = TRUE)
cmp2(t1, t1f)
cmp2(t1v, t1vf)

### and now maximally selected statistics
n <- 100
B <- 1000
x <- round(runif(n), 2)
y <- rnorm(n, mean = x < .5, sd = 2.6)
y2 <- runif(n)
blk <- gl(4, n/4)
ux <- sort(unique(x))
ix <- unclass(cut(x, breaks = c(-Inf, ux[-length(ux)] + diff(ux) / 2, Inf)))

(mt <- maxstat_test(y ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, B = 1000)
(tst <- doTest(lev, teststat = "maximum"))
ux[tst$index]
(tst <- doTest(lev, teststat = "quadratic"))
ux[tst$index]
lev <- LinStatExpCov(ix, y, B = 1000, varonly = TRUE)
(tst <- doTest(lev, teststat = "maximum"))
ux[tst$index]
(tst <- doTest(lev, teststat = "quadratic"))
ux[tst$index]

(mt <- maxstat_test(y ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, block = blk, B = 1000)
(tst <- doTest(lev, teststat = "maximum"))
ux[tst$index]
(tst <- doTest(lev, teststat = "quadratic"))
ux[tst$index]
lev <- LinStatExpCov(ix, y, block = blk, B = 1000, varonly = TRUE)
try(tst <- doTest(lev, teststat = "maximum"))

(mt <- maxstat_test(y + y2 ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000)
(tst <- doTest(lev, teststat = "maximum"))
ux[tst$index]
(tst <- doTest(lev, teststat = "quadratic"))
ux[tst$index]
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000, varonly = TRUE)
(tst <- doTest(lev, teststat = "maximum"))
ux[tst$index]
(tst <- doTest(lev, teststat = "quadratic"))
ux[tst$index]

(mt <- maxstat_test(y + y2 ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 1000)
(tst <- doTest(lev, teststat = "maximum"))
ux[tst$index]
(tst <- doTest(lev, teststat = "quadratic"))
ux[tst$index]
lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 1000, varonly = TRUE)
try(tst <- doTest(lev, teststat = "maximum"))

x <- sample(gl(5, n))
y <- rnorm(length(x), mean = x %in% levels(x)[c(1, 3, 5)], sd = 4.5)
y2 <- runif(length(x))
ix <- unclass(x)
blk <- gl(5, n)

(mt <- maxstat_test(y ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, B = 1000)
(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))
(tst <- doTest(lev, teststat = "quadratic", ordered = FALSE))
lev <- LinStatExpCov(ix, y, B = 1000, varonly = TRUE)
(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))
(tst <- doTest(lev, teststat = "quadratic", ordered = FALSE))

(mt <- maxstat_test(y ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, block = blk, B = 1000)
(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))
(tst <- doTest(lev, teststat = "quadratic", ordered = FALSE))
lev <- LinStatExpCov(ix, y, block = blk, B = 1000, varonly = TRUE)
try(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))

(mt <- maxstat_test(y + y2 ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000)
(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))
(tst <- doTest(lev, teststat = "quadratic", ordered = FALSE))
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000, varonly = TRUE)
(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))
(tst <- doTest(lev, teststat = "quadratic", ordered = FALSE))

(mt <- maxstat_test(y + y2 ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 50)
(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))
(tst <- doTest(lev, teststat = "quadratic", ordered = FALSE))

xx <- factor(x == levels(x)[tst$index == 1])
(it <- independence_test(y + y2 ~ xx | blk, teststat = "quadratic"))

lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 1000, varonly = TRUE)
try(tst <- doTest(lev, teststat = "maximum", ordered = FALSE))
