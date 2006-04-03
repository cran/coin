
### Regression tests for the K sample problem, i.e.,
### testing the independence of a numeric variable
### `y' and a factor `x' (possibly blocked)

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### generate data
dat <- data.frame(x = gl(4, 25), y = rnorm(100), block = gl(5, 20))[sample(1:100, 50),]

### Kruskal-Wallis Test

### asymptotic distribution
ptwo <- kruskal.test(y ~ x, data = dat)$p.value

stopifnot(isequal(pvalue(kruskal_test(y ~ x, data = dat)), ptwo))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
                                   teststat = "quad",
    ytrafo = function(data) trafo(data, numeric_trafo = rank))), ptwo))

### approximated distribution
rtwo <- pvalue(kruskal_test(y ~ x, data = dat, distribution = "approx")) / ptwo
stopifnot(all(rtwo > 0.9 & 
              rtwo < 1.1))

### <FIXME> add block examples </FIXME>

### sanity checks
try(kruskal_test(x ~ y, data = dat))
try(kruskal_test(x ~ y | y, data = dat))

### Fligner Test 

### asymptotic distribution
ptwo <- fligner.test(y ~ x, data = dat)$p.value                       

stopifnot(isequal(pvalue(fligner_test(y ~ x, data = dat)), ptwo))

dat$yy <- dat$y - tapply(dat$y, dat$x, median)[dat$x]
stopifnot(isequal(pvalue(oneway_test(yy ~ x, data = dat, distribution = "asympt", 
                                   teststat = "quad",
    ytrafo = function(data) trafo(data, numeric_trafo = fligner_trafo))), ptwo))

### approximated distribution
rtwo <- pvalue(fligner_test(y ~ x, data = dat, distribution = "approx")) / ptwo
stopifnot(all(rtwo > 0.9 &                     
              rtwo < 1.1))                     

### <FIXME> add block examples </FIXME>

### sanity checks
try(fligner_test(x ~ y, data = dat)) 
try(fligner_test(x ~ y | y, data = dat))

