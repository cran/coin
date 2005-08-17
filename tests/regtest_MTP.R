
### Regression tests for multiple adjustments

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### example from Westfall & Wolfinger (1997), Table 4
tab <- as.table(matrix(c(12, 1, 3, 8, 17, 12, 9, 9, 16, 24), nrow = 2,
    dimnames = list(group = c("Placebo", "Active"),
                    response = c("Very Poor", "Poor", "Fair", "Good",
                                 "Excellent"))))
df <- coin:::table2df(tab)

it <- independence_test(response ~ group, data = df,
                        distr = approximate(B = 100000))

### Table 5, last column: OK
pvalue(it, method = "step-down")

### example from Westfall & Wolfinger (1997), Table 2
df <- data.frame(group = factor(c(rep("Control", 50), rep("Treatment", 48))),
                 V1 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 5), rep(1, 5)),
                 V2 = c(rep(0, 50 - 4), rep(1, 4), rep(0, 48 - 3), rep(1, 3)),
                 V3 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 4), rep(1, 4)),
                 V4 = c(rep(0, 50 - 6), rep(1, 6), rep(0, 48 - 4), rep(1, 4)))

### alternative: less
it <- independence_test(V1 + V2 + V3 + V4 ~ group, data = df,
                        distr = approximate(B = 100000), alt = "less")

### page 4, 2nd column: adjusted p-value = 0.03665 for V1
pvalue(it, method = "discrete")

### alternative: less       
it <- independence_test(V1 + V2 + V3 + V4 ~ group, data = df,
                        distr = approximate(B = 100000), alt = "two")

### page 5, 1st column: adjusted p-value = 0.05261 for V1
pvalue(it, method = "discrete")

