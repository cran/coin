
### Regression tests for fixed bugs

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### I() returns objects of class AsIs which caused an error in `trafo'
df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = gl(2, 50))
independence_test(I(x1 / x2) ~ x3, data = df)
 