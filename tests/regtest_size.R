### Regression tests for actual test size

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)


### Berger and Ivanova (2001, p. 352, Tab 14.1)
tab <- as.table(matrix(c(12, 3, 10, 19), ncol = 2))

### Berger and Ivanova (2001, p. 355, Tab 14.2)
ct <- cmh_test(tab, distribution = "exact", alternative = "less",
               scores = list(Var2 = 0:1))
stopifnot(isequal(round(pvalue(ct), 4), 0.0049))
stopifnot(isequal(round(size(ct, alpha = 0.05), 4), 0.0273))

### Additional results
stopifnot(isequal(round(size(ct, alpha = 0.05, type = "mid-p-value"), 4), 0.0273))

ct <- cmh_test(tab, distribution = "approximate", alternative = "less",
               scores = list(Var2 = 0:1))
stopifnot(isequal(round(pvalue(ct), 4), 0.0045))
stopifnot(isequal(round(size(ct, alpha = 0.05), 4), 0.0254))
stopifnot(isequal(round(size(ct, alpha = 0.05, type = "mid-p-value"), 4), 0.0254))


### Berger and Ivanova (2001, p. 358, Tab 14.3)
tab <- as.table(matrix(c(12, 3, 3, 7, 7, 12), ncol = 3))

### Berger and Ivanova (2001, p. 360, Fig. 14.3, Fig 14.4)
ct_0 <- cmh_test(tab, distribution = "exact", alternative = "less",
                 scores = list(Var2 = c(0, 0, 1)))
stopifnot(isequal(round(pvalue(ct_0), 4), 0.1116))
stopifnot(isequal(round(size(ct_0, alpha = 0.05), 4), 0.0333))

ct_0.5 <- cmh_test(tab, distribution = "exact", alternative = "less",
                   scores = list(Var2 = c(0, 0.5, 1)))
stopifnot(isequal(round(pvalue(ct_0.5), 4), 0.0126))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05), 4), 0.0296))

ct_1 <- cmh_test(tab, distribution = "exact", alternative = "less",
                 scores = list(Var2 = c(0, 1, 1)))
stopifnot(isequal(round(pvalue(ct_1), 4), 0.0049))
stopifnot(isequal(round(size(ct_1, alpha = 0.05), 4), 0.0273))

### Additional results
stopifnot(isequal(round(size(ct_0, alpha = 0.05, type = "mid-p-value"), 4), 0.0333))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05, type = "mid-p-value"), 4), 0.0619))
stopifnot(isequal(round(size(ct_1, alpha = 0.05, type = "mid-p-value"), 4), 0.0273))

ct_0 <- cmh_test(tab, distribution = "approximate", alternative = "less",
                 scores = list(Var2 = c(0, 0, 1)))
stopifnot(isequal(round(pvalue(ct_0), 4), 0.1072))
stopifnot(isequal(round(size(ct_0, alpha = 0.05), 4), 0.0315))
stopifnot(isequal(round(size(ct_0, alpha = 0.05, type = "mid-p-value"), 4), 0.0315))

ct_0.5 <- cmh_test(tab, distribution = "approximate", alternative = "less",
                   scores = list(Var2 = c(0, 0.5, 1)))
stopifnot(isequal(round(pvalue(ct_0.5), 4), 0.0116))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05), 4), 0.0307))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05, type = "mid-p-value"), 4), 0.0614))

ct_1 <- cmh_test(tab, distribution = "approximate", alternative = "less",
                 scores = list(Var2 = c(0, 1, 1)))
stopifnot(isequal(round(pvalue(ct_1), 4), 0.0055))
stopifnot(isequal(round(size(ct_1, alpha = 0.05), 4), 0.0263))
stopifnot(isequal(round(size(ct_1, alpha = 0.05, type = "mid-p-value"), 4), 0.0263))


### Berger and Ivanova (2001, p. 364, Tab. 14.4)
tab <- as.table(array(c(4, 0, 14, 19, 11, 1, 15, 25), dim = c(2, 2, 2)))

### Berger and Ivanova (2001, p. 367, Fig. 14.8)
ct <- cmh_test(tab, distribution = "exact", alternative = "less",
               scores = list(Var2 = 0:1))
stopifnot(isequal(round(pvalue(ct),             4), 0.0001))
stopifnot(isequal(round(size(ct, alpha = 0.05), 4), 0.0229))

### Additional results
stopifnot(isequal(round(size(ct, alpha = 0.05, type = "mid-p-value"), 4), 0.0229))

ct <- cmh_test(tab, distribution = "approximate", alternative = "less",
               scores = list(Var2 = 0:1))
stopifnot(isequal(round(pvalue(ct), 4), 0.0001))
stopifnot(isequal(round(size(ct, alpha = 0.05), 4), 0.0229))
stopifnot(isequal(round(size(ct, alpha = 0.05, type = "mid-p-value"), 4), 0.0749))


### Berger and Ivanova (2002, p. 269)
tab <- as.table(matrix(c(11, 7, 2, 7, 2, 6), ncol = 3))

### Berger and Ivanova (2002, p. 277, Tab. 14.2, last line, first line)
ct_0 <- cmh_test(tab, distribution = "exact", alternative = "less",
                 scores = list(Var2 = c(0, 0, 1)))
stopifnot(isequal(round(pvalue(ct_0), 3), 0.228))
stopifnot(isequal(round(size(ct_0, alpha = 0.05), 3), 0.005))

ct_0.5 <- cmh_test(tab, distribution = "exact", alternative = "less",
                   scores = list(Var2 = c(0, 0.5, 1)))
stopifnot(isequal(round(pvalue(ct_0.5), 3), 0.038))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05), 3), 0.038))

ct_1 <- cmh_test(tab, distribution = "exact", alternative = "less",
                 scores = list(Var2 = c(0, 1, 1)))
stopifnot(isequal(round(pvalue(ct_1), 3), 0.028))
stopifnot(isequal(round(size(ct_1, alpha = 0.05), 3), 0.028))

### Additional results
stopifnot(isequal(round(size(ct_0, alpha = 0.05, type = "mid-p-value"), 3), 0.055))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05, type = "mid-p-value"), 3), 0.038))
stopifnot(isequal(round(size(ct_1, alpha = 0.05, type = "mid-p-value"), 3), 0.028))

ct_0   <- cmh_test(tab, distribution = "approximate", alternative = "less",
                   scores = list(Var2 = c(0, 0, 1)))
stopifnot(isequal(round(pvalue(ct_0), 3), 0.233))
stopifnot(isequal(round(size(ct_0, alpha = 0.05), 3), 0.006))
stopifnot(isequal(round(size(ct_0, alpha = 0.05, type = "mid-p-value"), 3), 0.058))

ct_0.5 <- cmh_test(tab, distribution = "approximate", alternative = "less",
                   scores = list(Var2 = c(0, 0.5, 1)))
stopifnot(isequal(round(pvalue(ct_0.5), 3), 0.039))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05), 3), 0.039))
stopifnot(isequal(round(size(ct_0.5, alpha = 0.05, type = "mid-p-value"), 3), 0.039))

ct_1   <- cmh_test(tab, distribution = "approximate", alternative = "less",
                   scores = list(Var2 = c(0, 1, 1)))
stopifnot(isequal(round(pvalue(ct_1), 3), 0.027))
stopifnot(isequal(round(size(ct_1, alpha = 0.05), 3), 0.027))
stopifnot(isequal(round(size(ct_1, alpha = 0.05, type = "mid-p-value"), 3), 0.027))


### Neuhäuser (2012, p. 30, Tab. 2.7)
y_1 <- rnorm(10); x_1 <- factor(rep(1:2, c(5, 5)))
y_2 <- rnorm(12); x_2 <- factor(rep(1:2, c(6, 6)))
y_3 <- rnorm(14); x_3 <- factor(rep(1:2, c(7, 7)))
y_4 <- rnorm(16); x_4 <- factor(rep(1:2, c(8, 8)))
y_5 <- rnorm(18); x_5 <- factor(rep(1:2, c(9, 9)))
y_6 <- rnorm(20); x_6 <- factor(rep(1:2, c(10, 10)))
y_7 <- rnorm(13); x_7 <- factor(rep(1:2, c(8, 5)))
y_8 <- rnorm(16); x_8 <- factor(rep(1:2, c(9, 7)))
y_9 <- rnorm(15); x_9 <- factor(rep(1:2, c(10, 5)))

wt_1 <- wilcox_test(y_1 ~ x_1, distribution = "exact")
stopifnot(isequal(round(size(wt_1, alpha = 0.05), 4), 0.0317))

wt_2 <- wilcox_test(y_2 ~ x_2, distribution = "exact")
stopifnot(isequal(round(size(wt_2, alpha = 0.05), 4), 0.0411))

wt_3 <- wilcox_test(y_3 ~ x_3, distribution = "exact")
stopifnot(isequal(round(size(wt_3, alpha = 0.05), 4), 0.0379))

wt_4 <- wilcox_test(y_4 ~ x_4, distribution = "exact")
stopifnot(isequal(round(size(wt_4, alpha = 0.05), 4), 0.0499))

wt_5 <- wilcox_test(y_5 ~ x_5, distribution = "exact")
stopifnot(isequal(round(size(wt_5, alpha = 0.05), 4), 0.0400))

wt_6 <- wilcox_test(y_6 ~ x_6, distribution = "exact")
stopifnot(isequal(round(size(wt_6, alpha = 0.05), 4), 0.0433))

wt_7 <- wilcox_test(y_7 ~ x_7, distribution = "exact")
stopifnot(isequal(round(size(wt_7, alpha = 0.05), 4), 0.0451))

wt_8 <- wilcox_test(y_8 ~ x_8, distribution = "exact")
stopifnot(isequal(round(size(wt_8, alpha = 0.05), 4), 0.0418))

wt_9 <- wilcox_test(y_9 ~ x_9, distribution = "exact")
stopifnot(isequal(round(size(wt_9, alpha = 0.05), 4), 0.0400))

### Additional results
stopifnot(isequal(round(size(wt_1, alpha = 0.05, type = "mid-p-value"), 4), 0.0556))
stopifnot(isequal(round(size(wt_2, alpha = 0.05, type = "mid-p-value"), 4), 0.0411))
stopifnot(isequal(round(size(wt_3, alpha = 0.05, type = "mid-p-value"), 4), 0.0530))
stopifnot(isequal(round(size(wt_4, alpha = 0.05, type = "mid-p-value"), 4), 0.0499))
stopifnot(isequal(round(size(wt_5, alpha = 0.05, type = "mid-p-value"), 4), 0.0503))
stopifnot(isequal(round(size(wt_6, alpha = 0.05, type = "mid-p-value"), 4), 0.0524))
stopifnot(isequal(round(size(wt_7, alpha = 0.05, type = "mid-p-value"), 4), 0.0451))
stopifnot(isequal(round(size(wt_8, alpha = 0.05, type = "mid-p-value"), 4), 0.0549))
stopifnot(isequal(round(size(wt_9, alpha = 0.05, type = "mid-p-value"), 4), 0.0553))

wt_1 <- wilcox_test(y_1 ~ x_1, distribution = "approximate")
stopifnot(isequal(round(size(wt_1, alpha = 0.05), 4), 0.0353))
stopifnot(isequal(round(size(wt_1, alpha = 0.05, type = "mid-p-value"), 4), 0.0582))

wt_2 <- wilcox_test(y_2 ~ x_2, distribution = "approximate")
stopifnot(isequal(round(size(wt_2, alpha = 0.05), 4), 0.0388))
stopifnot(isequal(round(size(wt_2, alpha = 0.05, type = "mid-p-value"), 4), 0.0388))

wt_3 <- wilcox_test(y_3 ~ x_3, distribution = "approximate")
stopifnot(isequal(round(size(wt_3, alpha = 0.05), 4), 0.0386))
stopifnot(isequal(round(size(wt_3, alpha = 0.05, type = "mid-p-value"), 4), 0.0530))

wt_4 <- wilcox_test(y_4 ~ x_4, distribution = "approximate")
stopifnot(isequal(round(size(wt_4, alpha = 0.05), 4), 0.0373))
stopifnot(isequal(round(size(wt_4, alpha = 0.05, type = "mid-p-value"), 4), 0.0508))

wt_5 <- wilcox_test(y_5 ~ x_5, distribution = "approximate")
stopifnot(isequal(round(size(wt_5, alpha = 0.05), 4), 0.0498))
stopifnot(isequal(round(size(wt_5, alpha = 0.05, type = "mid-p-value"), 4), 0.0498))

wt_6 <- wilcox_test(y_6 ~ x_6, distribution = "approximate")
stopifnot(isequal(round(size(wt_6, alpha = 0.05), 4), 0.0477))
stopifnot(isequal(round(size(wt_6, alpha = 0.05, type = "mid-p-value"), 4), 0.0477))

wt_7 <- wilcox_test(y_7 ~ x_7, distribution = "approximate")
stopifnot(isequal(round(size(wt_7, alpha = 0.05), 4), 0.0414))
stopifnot(isequal(round(size(wt_7, alpha = 0.05, type = "mid-p-value"), 4), 0.0414))

wt_8 <- wilcox_test(y_8 ~ x_8, distribution = "approximate")
stopifnot(isequal(round(size(wt_8, alpha = 0.05), 4), 0.0438))
stopifnot(isequal(round(size(wt_8, alpha = 0.05, type = "mid-p-value"), 4), 0.0438))

wt_9 <- wilcox_test(y_9 ~ x_9, distribution = "approximate")
stopifnot(isequal(round(size(wt_9, alpha = 0.05), 4), 0.0402))
stopifnot(isequal(round(size(wt_9, alpha = 0.05, type = "mid-p-value"), 4), 0.0549))
