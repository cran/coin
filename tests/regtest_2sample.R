
### Regression tests for the 2 sample problem, i.e.,
### testing the independence of a numeric variable
### `y' and a binary factor `x' (possibly blocked)

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### generate data
dat <- data.frame(x = gl(2, 50), y = rnorm(100), block = gl(5, 20))[sample(1:100,
75),]

### Wilcoxon Mann-Whitney Rank Sum Test

### asymptotic distribution
ptwo <- wilcox.test(y ~ x, data = dat, correct = FALSE, exact = FALSE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less", 
                     correct = FALSE, exact = FALSE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater", 
                        correct = FALSE, exact = FALSE)$p.value

stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "less")),
                  pless))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater")), 
                  pgreater))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank))), ptwo))
### check direct supply of a function via ytrafo
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = rank)), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "less")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "greater")), pgreater))

### exact distribution
ptwo <- wilcox.test(y ~ x, data = dat, exact = TRUE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less", exact = TRUE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater", 
                        exact = TRUE)$p.value

stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, distribution = "exact")), 
                  ptwo))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "less", 
                             distribution = "exact")), pless))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater",
                             distribution = "exact")), pgreater))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank))), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "less")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "greater")), pgreater))

### approximated distribution
rtwo <- pvalue(wilcox_test(y ~ x, data = dat, distribution = "approx")) / ptwo
rless <- pvalue(wilcox_test(y ~ x, data = dat, alternative = "less",
                   distribution = "approx")) / pless
rgreater <- pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater",
                   distribution = "approx")) / pgreater
stopifnot(all(c(rtwo, rless, rgreater) > 0.9 & 
              c(rtwo, rless, rgreater) < 1.1))

### <FIXME> add block examples </FIXME>

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact"))

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt", 
                   alternative = "less"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx", 
                   alternative = "less"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact", 
                   alternative = "less"))

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt", 
                   alternative = "greater"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx", 
                   alternative = "greater"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact", 
                   alternative = "greater"))

### sanity checks
try(wilcox_test(x ~ y, data = dat))
try(wilcox_test(x ~ y | y, data = dat))

### Ansari-Bradley Test 

### asymptotic distribution
ptwo <- ansari.test(y ~ x, data = dat, correct = FALSE, exact = FALSE)$p.value
pless <- ansari.test(y ~ x, data = dat, alternative = "less", 
                     correct = FALSE, exact = FALSE)$p.value
pgreater <- ansari.test(y ~ x, data = dat, alternative = "greater", 
                        correct = FALSE, exact = FALSE)$p.value

stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "less")),
                  pless))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "greater")), 
                  pgreater))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo))), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "greater")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "less")), pgreater))

### exact distribution
ptwo <- ansari.test(y ~ x, data = dat, exact = TRUE)$p.value
pless <- ansari.test(y ~ x, data = dat, alternative = "less", exact = TRUE)$p.value
pgreater <- ansari.test(y ~ x, data = dat, alternative = "greater", 
                        exact = TRUE)$p.value

### <FIXME>: Definition of two-sided P-values! </FIXME>
(isequal(pvalue(ansari_test(y ~ x, data = dat, distribution = "exact")), 
                  ptwo))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "less", 
                             distribution = "exact")), pless))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "greater",
                             distribution = "exact")), pgreater))

### <FIXME>: Definition of two-sided P-values! </FIXME>
(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo))), ptwo))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "greater")), pless))
stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "less")), pgreater))

### approximated distribution
rtwo <- pvalue(ansari_test(y ~ x, data = dat, distribution = "approx")) / ptwo
rless <- pvalue(ansari_test(y ~ x, data = dat, alternative = "less",
                   distribution = "approx")) / pless
rgreater <- pvalue(ansari_test(y ~ x, data = dat, alternative = "greater",
                   distribution = "approx")) / pgreater
### <FIXME> ??? </FIXME>
(all(c(rtwo, rless, rgreater) > 0.9 & 
              c(rtwo, rless, rgreater) < 1.1))

### <FIXME> add block examples </FIXME>

### sanity checks
try(ansari_test(x ~ y, data = dat))
try(ansari_test(x ~ y | y, data = dat))

### the remaining three candidates
### asymptotic distribution
ptwo <- wilcox.test(y ~ x, data = dat, correct = FALSE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less",
                     correct = FALSE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater", 
                        correct = FALSE)$p.value

oneway_test(y ~ x, dat = dat)
oneway_test(y ~ x, dat = dat, alternative = "less")
oneway_test(y ~ x, dat = dat, alternative = "greater")

normal_test(y ~ x, dat = dat)
normal_test(y ~ x, dat = dat, alternative = "less")
normal_test(y ~ x, dat = dat, alternative = "greater")

median_test(y ~ x, dat = dat)
median_test(y ~ x, dat = dat, alternative = "less")
median_test(y ~ x, dat = dat, alternative = "greater")

### confidence intervals, cf Bauer 1972
location <- data.frame(y = c(6, 20, 27, 38, 46, 51, 54, 57,
                             10, 12, 15, 21, 32, 40, 41, 45),
                       x = gl(2, 8))

ci <- confint(normal_test(y ~ x, data = location, 
                          conf.int = TRUE, di = "ex"))
stopifnot(isequal(ci$conf.int, c(-6, 30)))
stopifnot(isequal(ci$estimate, 11))

wt <- wilcox.test(y ~ x, data = location, 
                  conf.int = TRUE)
ci <- confint(wilcox_test(y ~ x, data = location, 
                          conf.int = TRUE, di = "ex"))
stopifnot(isequal(wt$confint, ci$confint))
stopifnot(isequal(wt$estimate, ci$estimate))

scale <- data.frame(y = c(-101, -35, -13, 10, 130, 236, 370, 556,
                          -145, -140, -40, -30, 2, 27, 68, 290),
                    x = gl(2, 8))

ci <- confint(ansari_test(y ~ x, data = scale, conf.int = TRUE, 
                          di = "ex", conf.level = 0.988))
stopifnot(isequal(ci$conf.int, c(10, 556) / c(68, 27)))
stopifnot(isequal(ci$estimate, mean(c(35/30, 370 / 290))))

### ties handling
y1 <- c(14 , 18 , 2 , 4 , -5 , 14 , -3 , -1 , 1 , 6 , 3 , 3)
x1 <- c(8 , 26 , -7 , -1 , 2 , 9 , 0 , -4 , 13 , 3 , 3 , 4)
pvalue(wilcoxsign_test(y1~x1,alter="greater",dist=exact(), 
                       zero.method = "Wilcoxon"))
pvalue(wilcoxsign_test(y1~x1,alter="greater",dist=exact()))
