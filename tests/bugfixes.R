
### Regression tests for fixed bugs

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### I() returns objects of class AsIs which caused an error in `trafo'
df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = gl(2, 50))
independence_test(I(x1 / x2) ~ x3, data = df)

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


