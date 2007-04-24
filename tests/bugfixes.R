
### Regression tests for fixed bugs

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### check if doxygen documentation is there
stopifnot(nchar(system.file("documentation", "html", "index.html", 
                            package = "coin")) > 29) 

### I() returns objects of class AsIs which caused an error in `trafo'
df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = gl(2, 50))
independence_test(I(x1 / x2) ~ x3, data = df)
independence_test(I(x1 < 0) ~ x3, data = df)

### expectation was wrong when varonly = TRUE in case both
### xtrafo and ytrafo were multivariate 
if (require(multcomp)) {
    df <- data.frame(x = runif(30), y = runif(30), z = gl(3, 10))
    a <- independence_test(x + y ~ z, data = df,
         distribution = approximate(B = 19999),
         xtrafo = function(data) trafo(data, factor_trafo = function(x)
             model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))))
    b <- independence_test(x + y ~ z, data = df,
         xtrafo = function(data) trafo(data, factor_trafo = function(x)
             model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))))
    isequal(expectation(a), expectation(b))
}


### `statistic' for linear and standardized statistics was wrong in case of 
### scores
data("jobsatisfaction")
stopifnot(unique(dim(statistic(lbl_test(jobsatisfaction), "linear"))) == 1)
stopifnot(unique(dim(statistic(lbl_test(jobsatisfaction), "standardized"))) == 1)


### support() failed in most cases
df <- data.frame(x = runif(20), y = runif(20), z = gl(2, 10))
support(independence_test(x ~ z, data = df))
support(independence_test(x ~ z, data = df, teststat = "quad"))
ite <- independence_test(I(round(x, 1)) ~ z, data = df, dist = exact())
ae <- support(ite)
de <- sapply(ae, function(x) dperm(ite, x))
sum(de)
ita <- independence_test(I(round(x, 1)) ~ z, data = df, 
                         dist = approximate(B = 100000))
aa <- support(ita)
da <- sapply(aa, function(x) dperm(ita, x))
sum(da)
mean(round(ae, 10) %in% aa)

plot(aa, da, type = "s", lty = 1)
lines(ae, de, type = "s", lty = 2)
itas <- independence_test(I(round(x, 1)) ~ z, data = df)
lines(ae[-1], diff(sapply(ae, function(x) pperm(itas, x))), lty = 3)
legend("topleft", lty = 1:3, legend = c("approx", "exact", "asympt"), bty = "n")

### check correct handling of multiple censoring indicators (in modeltools)
### was never wrong, just in case...
data("photocar", package = "coin")
i1 <- independence_test(Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group,
                  data = photocar)
i2 <- independence_test(Surv(time, event) ~ group, data = photocar)
i3 <- independence_test(Surv(dmin, tumor) ~ group, data = photocar)

stopifnot(max(abs(statistic(i1, "standardized")[,1] - 
                  statistic(i2, "stand"))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(statistic(i1, "standardized")[,2] - 
                  statistic(i3, "stand"))) < sqrt(.Machine$double.eps))

### check new var_trafo argument
x <- rnorm(20)
y <- gl(2, 10)
a <- trafo(data.frame(x = x, y = y), numeric_trafo = normal_trafo)
b <- trafo(data.frame(x = x, y = y), var_trafo = list(x = normal_trafo))
stopifnot(isequal(a, b))

### check for multiple ordered factors
mydf <- data.frame(x = ordered(gl(4, 5)), y = ordered(gl(5, 4)), 
                   z = rnorm(20))
it1 <- independence_test(x + y ~ z , data = mydf)
stopifnot(isequal(drop(statistic(it1, "linear")), 
          c(statistic(independence_test(x ~ z , data = mydf), "linear"),
            statistic(independence_test(y ~ z , data = mydf), "linear"))))
it1 <- independence_test(x + z ~ y , data = mydf)
stopifnot(isequal(drop(statistic(it1, "linear")),
          c(statistic(independence_test(x ~ y , data = mydf), "linear"),
            statistic(independence_test(z ~ y , data = mydf), "linear"))))
it1 <- independence_test(z ~ x + y , data = mydf)
stopifnot(isequal(drop(statistic(it1, "linear")),
          c(statistic(independence_test(z ~ x , data = mydf), "linear"),
            statistic(independence_test(z ~ y , data = mydf), "linear"))))

### NA's and weights
mydf <- data.frame(x = 1:10, y = gl(2, 5), w = rep(2, 10))
s <- statistic(independence_test(x ~ y, data = mydf, weights = ~ w), "linear")
stopifnot(s == 30)
mydf$x[1] <- NA
s <- statistic(independence_test(x ~ y, data = mydf, weights = ~ w), "linear")
stopifnot(s == 28)

### two observations only
mydf <- data.frame(x = 1:10, y = factor(rep(c(1, 2), 5)))
independence_test(y ~ x, data = mydf, subset = c(1, 6))
independence_test(y ~ x, data = mydf, subset = c(1, 2))
try(independence_test(y ~ x, data = mydf, subset = 1))

### names of expectation and covariance
YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
                             42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
                             38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
                             31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                  site = factor(c(rep("I", 10), rep("II", 10),
                                  rep("III", 10), rep("IV", 10))))

it <- independence_test(length ~ site, data = YOY,
    ytrafo = function(data) trafo(data, numeric_trafo = rank),
    teststat = "quad")
expectation(it)
covariance(it)

mydf <- data.frame(x = rnorm(10), y = rnorm(10), z = gl(2, 5))
it <- independence_test(x + y ~ z, data = mydf)
statistic(it, "linear")
expectation(it)
covariance(it)

### maxstat_trafo
n <- seq(from = 5, to = 100, by = 1)
for (i in n) {
   x <- round(rnorm(i) * 10, 1)
   xm <- maxstat_trafo(x)
   stopifnot(min(c(mean(xm[,1]), 1 - mean(xm[,ncol(xm)])) - 0.1) >
             -.Machine$double.eps)
}

### formula evaluation in `parent.frame()', spotted by Z
foo <- function(x, y) independence_test(y ~ x)
a <- 1:10
b <- 1:10
foo(a, b)
x <- 1
y <- 1
foo(a, b)

### factors with only one level
dat <- data.frame(y = rnorm(100), x1 = runif(100), x2 = factor(rep(0, 100)))
try(independence_test(y ~ x1  + x2, data = dat))

### user specified g: names, MC
me <- as.table(matrix(c( 6,  8, 10,
               32, 47, 20), byrow = TRUE, nrow = 2,
    dimnames = list(group = c("In situ", "Control"),
                    genotype = c("AA", "AG", "GG"))))
medf <- as.data.frame(me)

add <- c(0, 1, 2)
dom <- c(0, 1, 1)
rez <- c(0, 0, 1)
g <- function(x) {
    x <- unlist(x)
    cbind(add[x], dom[x], rez[x])
}
it <- independence_test(group ~ genotype, 
    data = medf, weights = ~ Freq, xtrafo = g)
statistic(it, "linear")

it <- independence_test(group ~ genotype,
    data = medf, weights = ~ Freq, xtrafo = g,
    distribution = approximate(B = 49999))
pvalue(it)

stopifnot(all.equal(statistic(independence_test(t(me), xtrafo = g), "linear"),
                    statistic(it, "linear")))
