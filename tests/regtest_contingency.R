
### Regression tests for the r x c x K problem, i.e.,
### testing the independence of a factor
### `y' and a factor factor `x' (possibly blocked)

set.seed(290875)
library(coin)
isequal <- coin:::isequal

thisversion <- paste(R.version$major, R.version$minor, sep = ".")

### generate data: 2 x 2 x K
dat <- data.frame(x = gl(2, 50), y = gl(2, 50)[sample(1:100)], 
                  block = gl(10, 10)[sample(1:100)])[sample(1:100, 75),]

### Pearsons Chisq Test, asymptotic distribution
ptwo <- chisq.test(table(dat$x, dat$y), correct = FALSE)$p.value

stopifnot(isequal(pvalue(chisq_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(chisq_test(table(dat$y, dat$x))), ptwo))

### Cochran-Mantel-Haenzel Test, asymptotic distribution
ptwo <- drop(mantelhaen.test(table(dat$x, dat$y, dat$block), 
                             correct = FALSE)$p.value)

stopifnot(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
stopifnot(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))


### generate data: r x c x K
dat <- data.frame(x = gl(4, 25), y = gl(4, 25)[sample(1:100)], 
                  block = gl(2, 50)[sample(1:100)])

### Cochran-Mantel-Haenzel Test, asymptotic distribution
### _is wrong_ in R < 2.1.0!!!
ptwo <- drop(mantelhaen.test(table(dat$y, dat$x, dat$block), 
                             correct = FALSE)$p.value)

if (compareVersion(thisversion, "2.1.0") >= 0) {
    stopifnot(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
    stopifnot(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))
}

### generate data: r x c x K
dat <- data.frame(x = gl(4, 25), y = gl(5, 20)[sample(1:100)], 
                  block = gl(2, 50)[sample(1:100)])

### Cochran-Mantel-Haenzel Test, asymptotic distribution
### _is wrong_!!!
ptwo <- drop(mantelhaen.test(table(dat$y, dat$x, dat$block),
                             correct = FALSE)$p.value)

if (compareVersion(thisversion, "2.1.0") >= 0) {
    stopifnot(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
    stopifnot(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))
}

### 2x2 table and maxstat
x <- c(rep(1,51), rep(2,49))
y <- factor(c(rep(0,49), rep(1,51)))[sample(1:100)]
stopifnot(isequal(as.vector(statistic(independence_test(table(x, y)))),
as.vector(statistic(maxstat_test(y ~ x )))))

### marginal homogeneity
rating <- c("low", "moderate", "high")
x <- as.table(matrix(c(20, 10,  5,
                       3, 30, 15, 
                       0,  5, 40),
                     ncol = 3, byrow = TRUE,
                     dimnames = list(Rater1 = rating, Rater2 = rating)))
### test statistic W_0 = 13.76
### see http://ourworld.compuserve.com/homepages/jsuebersax/mcnemar.htm
stopifnot(all.equal(round(statistic(mh_test(x)), 2), 13.76))


### see `comparison.R' for more regression tests