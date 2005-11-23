
### Regression tests for fixed bugs

set.seed(290875)
library(coin)
isequal <- coin:::isequal

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

