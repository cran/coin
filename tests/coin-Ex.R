### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign(".CheckExEnv", as.environment(2), pos = length(search())) # base
## add some hooks to label plot pages for base and grid graphics
setHook("plot.new", ".newplot.hook")
setHook("persp", ".newplot.hook")
setHook("grid.newpage", ".gridplot.hook")

assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       env = .CheckExEnv)
assign("..nameEx", "__{must remake R-ex/*.R}__", env = .CheckExEnv) # for now
assign("ptime", proc.time(), env = .CheckExEnv)
grDevices::postscript("coin-Examples.ps")
assign("par.postscript", graphics::par(no.readonly = TRUE), env = .CheckExEnv)
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
library('coin')

assign(".oldSearch", search(), env = .CheckExEnv)
assign(".oldNS", loadedNamespaces(), env = .CheckExEnv)
cleanEx(); ..nameEx <- "ContingencyTests"

### * ContingencyTests

flush(stderr()); flush(stdout())

### Name: ContingencyTests
### Title: Independence in General I x K x J Contingency Tables
### Aliases: chisq_test chisq_test.formula chisq_test.table
###   chisq_test.IndependenceProblem cmh_test.formula cmh_test.table
###   cmh_test.IndependenceProblem cmh_test lbl_test.formula lbl_test.table
###   lbl_test.IndependenceProblem lbl_test
### Keywords: htest

### ** Examples


## Don't show: 
    set.seed(290875)
## End Don't show

data(jobsatisfaction, package = "coin")

### for females only
chisq_test(as.table(jobsatisfaction[,,"Female"]), 
    distribution = approximate(B = 9999))

### both Income and Job.Satisfaction unordered
cmh_test(jobsatisfaction)

### both Income and Job.Satisfaction ordered, default scores
lbl_test(jobsatisfaction)

### both Income and Job.Satisfaction ordered, alternative scores
lbl_test(jobsatisfaction, scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                                        Income = c(3, 10, 20, 35)))

### the same, null distribution approximated
cmh_test(jobsatisfaction, scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                                        Income = c(3, 10, 20, 35)),
         distribution = approximate(B = 10000))




cleanEx(); ..nameEx <- "IndependenceTest"

### * IndependenceTest

flush(stderr()); flush(stdout())

### Name: IndependenceTest
### Title: General Independence Tests
### Aliases: independence_test independence_test.formula
###   independence_test.IndependenceProblem independence_test.table
### Keywords: htest

### ** Examples


data(asat, package = "coin")

### independence of asat and group via normal scores test
independence_test(asat ~ group, data = asat,

    ### exact null distribution
    distribution = "exact", 

    ### one-sided test
    alternative = "greater",

    ### apply normal scores to asat$asat
    ytrafo = function(data) trafo(data, numeric_trafo = normal_trafo),

    ### indicator matrix of 1st level of group
    xtrafo = function(data) trafo(data, factor_trafo = function(x)
        matrix(x == levels(x)[1], ncol = 1))
)

### same as
normal_test(asat ~ group, data = asat, distribution = "exact", 
            alternative = "greater")




cleanEx(); ..nameEx <- "LocationTests"

### * LocationTests

flush(stderr()); flush(stdout())

### Name: LocationTests
### Title: Independent Two- and K-Sample Location Tests
### Aliases: oneway_test oneway_test.formula
###   oneway_test.IndependenceProblem wilcox_test.formula
###   wilcox_test.IndependenceProblem wilcox_test normal_test.formula
###   normal_test.IndependenceProblem normal_test median_test.formula
###   median_test.IndependenceProblem median_test kruskal_test.formula
###   kruskal_test.IndependenceProblem kruskal_test
### Keywords: htest

### ** Examples


### Tritiated Water Diffusion Across Human Chorioamnion
### Hollander & Wolfe (1999), Table 4.1, page 110
water_transfer <- data.frame(
    pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
           1.15, 0.88, 0.90, 0.74, 1.21),
    age = factor(c(rep("At term", 10), rep("12-26 Weeks", 5))))

### Wilcoxon-Mann-Whitney test, cf. Hollander & Wolfe (1999), page 111
### exact p-value and confidence interval for the difference in location
### (At term - 12-26 Weeks)
wt <- wilcox_test(pd ~ age, data = water_transfer, 
            distribution = "exact", conf.int = TRUE)
print(wt)

### extract observed Wilcoxon statistic, i.e, the sum of the
### ranks for age = "12-26 Weeks"
statistic(wt, "linear")

### its expectation
expectation(wt)

### and variance
covariance(wt)

### and, finally, the exact two-sided p-value
pvalue(wt)

### Confidence interval for difference (12-26 Weeks - At term)
wilcox_test(pd ~ age, data = water_transfer, 
    xtrafo = function(data) 
        trafo(data, factor_trafo = function(x) as.numeric(x == levels(x)[2])),
        distribution = "exact", conf.int = TRUE)

### Permutation test, asymptotic p-value
oneway_test(pd ~ age, data = water_transfer)

### approximate p-value (with 99% confidence interval)
pvalue(oneway_test(pd ~ age, data = water_transfer, 
                 distribution = approximate(B = 9999)))
### exact p-value
pt <- oneway_test(pd ~ age, data = water_transfer, distribution = "exact")
pvalue(pt)

### plot density and distribution of the standardised 
### test statistic
layout(matrix(1:2, nrow = 2))
s <- support(pt)
d <- sapply(s, function(x) dperm(pt, x))
p <- sapply(s, function(x) pperm(pt, x))
plot(s, d, type = "S", xlab = "Teststatistic", ylab = "Density")
plot(s, p, type = "S", xlab = "Teststatistic", ylab = "Cumm. Probability")

### Length of YOY Gizzard Shad from Kokosing Lake, Ohio,
### sampled in Summer 1984, Hollander & Wolfe, Table 6.3, page 200
YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44, 
                             42, 60, 32, 42, 45, 58, 27, 51, 42, 52, 
                             38, 33, 26, 25, 28, 28, 26, 27, 27, 27, 
                             31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                  site = factor(c(rep("I", 10), rep("II", 10),
                                  rep("III", 10), rep("IV", 10))))

### Kruskal-Wallis test, approximate exact p-value
kw <- kruskal_test(length ~ site, data = YOY, 
             distribution = approximate(B = 9999))
kw
pvalue(kw)

### Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
### Hollander & Wolfe, 1999, page 244 
### (where Steel-Dwass results are given)
if (require(multcomp)) {

    NDWD <- oneway_test(length ~ site, data = YOY,
        ytrafo = function(data) trafo(data, numeric_trafo = rank),
        xtrafo = function(data) trafo(data, factor_trafo = function(x)
            model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
        teststat = "maxtype", distribution = approximate(B = 90000))

    ### global p-value
    print(pvalue(NDWD))

    ### sites (I = II) != (III = IV) at alpha = 0.01 (page 244)
    print(pvalue(NDWD, adjusted = TRUE))
}




cleanEx(); ..nameEx <- "MarginalHomogenityTest"

### * MarginalHomogenityTest

flush(stderr()); flush(stdout())

### Name: MarginalHomogenityTest
### Title: Marginal Homogenity Test
### Aliases: mh_test mh_test.table mh_test.formula mh_test.SymmetryProblem
### Keywords: htest

### ** Examples


### Opinions on Pre- and Extramarital Sex, Agresti (2002), page 421
opinions <- c("always wrong", "almost always wrong", 
              "wrong only sometimes", "not wrong at all")

PreExSex <- as.table(matrix(c(144, 33, 84, 126, 
                                2,  4, 14,  29, 
                                0,  2,  6,  25, 
                                0,  0,  1,  5), nrow = 4, 
                            dimnames = list(PremaritalSex = opinions,
                                            ExtramaritalSex = opinions)))

### treating response as nominal
mh_test(PreExSex)

### and as ordinal
mh_test(PreExSex, scores = list(response = 1:length(opinions)))




cleanEx(); ..nameEx <- "MaxstatTest"

### * MaxstatTest

flush(stderr()); flush(stdout())

### Name: MaxstatTest
### Title: Maximally Selected Statistics
### Aliases: maxstat_test maxstat_test.formula
###   maxstat_test.IndependenceProblem
### Keywords: htest

### ** Examples


data(treepipit, package = "coin")

maxstat_test(counts ~ coverstorey, data = treepipit)




cleanEx(); ..nameEx <- "PermutationDistribution"

### * PermutationDistribution

flush(stderr()); flush(stdout())

### Name: PermutationDistribution
### Title: Permutation Distribution of Conditional Independence Tests
### Aliases: pperm pperm-methods pperm,NullDistribution-method
###   pperm,IndependenceTest-method pperm,ScalarIndependenceTest-method
###   pperm,MaxTypeIndependenceTest-method
###   pperm,QuadTypeIndependenceTest-method qperm qperm-methods
###   qperm,NullDistribution-method qperm,IndependenceTest-method
###   qperm,ScalarIndependenceTest-method
###   qperm,MaxTypeIndependenceTest-method
###   qperm,QuadTypeIndependenceTest-method dperm dperm-methods
###   dperm,NullDistribution-method dperm,IndependenceTest-method
###   dperm,ScalarIndependenceTest-method
###   dperm,MaxTypeIndependenceTest-method
###   dperm,QuadTypeIndependenceTest-method support support-methods
###   support,NullDistribution-method support,IndependenceTest-method
###   support,ScalarIndependenceTest-method
###   support,MaxTypeIndependenceTest-method
###   support,QuadTypeIndependenceTest-method
### Keywords: methods htest

### ** Examples


### artificial 2-sample problem
df <- data.frame(y = rnorm(20), x = gl(2, 10))

### Ansari-Bradley test
at <- ansari_test(y ~ x, data = df, distribution = "exact")

### density of the exact distribution of the Ansari-Bradley statistic
dens <- sapply(support(at), dperm, object = at)

### 95% quantile
qperm(at, 0.95)

### one-sided p-value
pperm(at, statistic(at))




cleanEx(); ..nameEx <- "ScaleTests"

### * ScaleTests

flush(stderr()); flush(stdout())

### Name: ScaleTests
### Title: Independent Two- and K-Sample Scale Tests
### Aliases: ansari_test ansari_test.formula
###   ansari_test.IndependenceProblem fligner_test.formula
###   fligner_test.IndependenceProblem fligner_test
### Keywords: htest

### ** Examples


### Serum Iron Determination Using Hyland Control Sera
### Hollander & Wolfe (1999), page 147
sid <- data.frame(
    serum = c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
              101, 96, 97, 102, 107, 113, 116, 113, 110, 98,
              107, 108, 106, 98, 105, 103, 110, 105, 104,
              100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99),
    method = factor(gl(2, 20), labels = c("Ramsay", "Jung-Parekh")))

### Ansari-Bradley test, asymptotical p-value
ansari_test(serum ~ method, data = sid)

### exact p-value
ansari_test(serum ~ method, data = sid, distribution = "exact")

### Platelet Counts of Newborn Infants
### Hollander & Wolfe, Table 5.4, page 171
platalet_counts <- data.frame(
    counts = c(120, 124, 215, 90, 67, 95, 190, 180, 135, 399, 
               12, 20, 112, 32, 60, 40),
    treatment = factor(c(rep("Prednisone", 10), rep("Control", 6))))

### Lepage test, Hollander & Wolfe, page 172 
lt <- independence_test(counts ~ treatment, data = platalet_counts,
    ytrafo = function(data) trafo(data, numeric_trafo = function(x)       
        cbind(rank(x), ansari_trafo(x))),
    teststat = "quadtype", distribution = approximate(B = 9999))

lt

### where did the rejection come from? Use maximum statistic
### instead of a quadratic form
ltmax <- independence_test(counts ~ treatment, data = platalet_counts,
    ytrafo = function(data) trafo(data, numeric_trafo = function(x) 
        matrix(c(rank(x), ansari_trafo(x)), ncol = 2,
               dimnames = list(1:length(x), c("Location", "Scale")))),
    teststat = "maxtype")

### points to a difference in location
pvalue(ltmax, adjusted = TRUE)

### Funny: We could have used a simple Bonferroni procedure
### since the correlation between the Wilcoxon and Ansari-Bradley 
### test statistics is zero
covariance(ltmax)




cleanEx(); ..nameEx <- "SpearmanTest"

### * SpearmanTest

flush(stderr()); flush(stdout())

### Name: SpearmanTest
### Title: Spearman's Test on Independence
### Aliases: spearman_test spearman_test.formula
###   spearman_test.IndependenceProblem
### Keywords: htest

### ** Examples


spearman_test(CONT ~ INTG, data = USJudgeRatings)




cleanEx(); ..nameEx <- "SurvTest"

### * SurvTest

flush(stderr()); flush(stdout())

### Name: SurvTest
### Title: Independent Two- and K-Sample Tests for Censored Data
### Aliases: surv_test surv_test.formula surv_test.IndependenceProblem
### Keywords: htest

### ** Examples


### asymptotic tests for carcinoma data
data(ocarcinoma, package = "coin")
surv_test(Surv(time, event) ~ stadium, data = ocarcinoma)
survdiff(Surv(time, event) ~ stadium, data = ocarcinoma) 

### example data given in Callaert (2003)
exdata <- data.frame(time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
                     event = rep(TRUE, 15),
                     group = factor(c(rep(0, 7), rep(1, 8))))
### p = 0.0523
survdiff(Surv(time, event) ~ group, data = exdata)
### p = 0.0505
surv_test(Surv(time, event) ~ group, data = exdata, 
          distribution = exact())


cleanEx(); ..nameEx <- "SymmetryTests"

### * SymmetryTests

flush(stderr()); flush(stdout())

### Name: SymmetryTests
### Title: Symmetry Tests
### Aliases: friedman_test friedman_test.formula
###   friedman_test.SymmetryProblem wilcoxsign_test.formula
###   wilcoxsign_test.IndependenceProblem wilcoxsign_test
### Keywords: htest

### ** Examples


### Hollander & Wolfe (1999), Table 7.1, page 274
### Comparison of three methods ("round out", "narrow angle", and
###  "wide angle") for rounding first base. 
RoundingTimes <- data.frame(
    times = c(5.40, 5.50, 5.55,
              5.85, 5.70, 5.75,
              5.20, 5.60, 5.50,
              5.55, 5.50, 5.40,
              5.90, 5.85, 5.70,
              5.45, 5.55, 5.60,
              5.40, 5.40, 5.35,
              5.45, 5.50, 5.35,
              5.25, 5.15, 5.00,
              5.85, 5.80, 5.70,
              5.25, 5.20, 5.10,
              5.65, 5.55, 5.45,
              5.60, 5.35, 5.45,
              5.05, 5.00, 4.95,
              5.50, 5.50, 5.40,
              5.45, 5.55, 5.50,
              5.55, 5.55, 5.35,
              5.45, 5.50, 5.55,
              5.50, 5.45, 5.25,
              5.65, 5.60, 5.40,
              5.70, 5.65, 5.55,
              6.30, 6.30, 6.25),
    methods = factor(rep(c("Round Out", "Narrow Angle", "Wide Angle"), 22)),
    block = factor(rep(1:22, rep(3, 22))))

### classical global test
friedman_test(times ~ methods | block, data = RoundingTimes)

### parallel coordinates plot
matplot(t(matrix(RoundingTimes$times, ncol = 3, byrow = TRUE)), 
        type = "l", col = 1, lty = 1, axes = FALSE, ylab = "Time", 
        xlim = c(0.5, 3.5))
axis(1, at = 1:3, labels = levels(RoundingTimes$methods))
axis(2)

### where do the differences come from?
### Wilcoxon-Nemenyi-McDonald-Thompson test
### Hollander & Wolfe, page 295
if (require(multcomp)) {

    ### all pairwise comparisons
    rtt <- symmetry_test(times ~ methods | block, data = RoundingTimes,
         teststat = "maxtype",
         xtrafo = function(data)
             trafo(data, factor_trafo = function(x)
                 model.matrix(~ x - 1) 
             ),
         ytrafo = function(data)
             trafo(data, numeric_trafo = rank, block = RoundingTimes$block)
    )

    ### a global test, again
    print(pvalue(rtt))

    ### simultaneous P-values for all pair comparisons
    ### Wide Angle vs. Round Out differ (Hollander and Wolfe, 1999, page 296)
    print(pvalue(rtt, adjusted = TRUE))
}

### Strength Index of Cotton, Hollander & Wolfe (1999), Table 7.5, page 286
sc <- data.frame(block = factor(c(rep(1, 5), rep(2, 5), rep(3, 5))),
                 potash = ordered(rep(c(144, 108, 72, 54, 36), 3)),
                 strength = c(7.46, 7.17, 7.76, 8.14, 7.63,
                              7.68, 7.57, 7.73, 8.15, 8.00,
                              7.21, 7.80, 7.74, 7.87, 7.93))

### Page test
ft <- friedman_test(strength ~ potash | block, data = sc)
ft

### one-sided p-value
1 - pnorm(sqrt(statistic(ft)))

### approximate null distribution via Monte-Carlo
pvalue(friedman_test(strength ~ potash | block, data = sc, 
                     distribution = approximate(B = 9999)))




cleanEx(); ..nameEx <- "Transformations"

### * Transformations

flush(stderr()); flush(stdout())

### Name: Transformations
### Title: Functions for Data Transformations and Score Computations
### Aliases: trafo id_trafo ansari_trafo fligner_trafo normal_trafo
###   median_trafo consal_trafo maxstat_trafo logrank_trafo f_trafo
### Keywords: manip

### ** Examples


### dummy matrices, 2-sample problem (only one column)
f_trafo(y <- gl(2, 5))

### K-sample problem (K columns)
f_trafo(y <- gl(5, 2))

### normal scores
normal_trafo(x <- rnorm(10))

### and now together
trafo(data.frame(x = x, y = y), numeric_trafo = normal_trafo)

### maximally selected statistics
maxstat_trafo(rnorm(10))

### apply transformation blockwise (e.g. for Friedman test)
trafo(data.frame(y = 1:20), numeric_trafo = rank, block = gl(4, 5))




cleanEx(); ..nameEx <- "asat"

### * asat

flush(stderr()); flush(stdout())

### Name: asat
### Title: Toxicological Study on Female Wistar Rats
### Aliases: asat
### Keywords: datasets

### ** Examples


data(asat, package = "coin")

### proof-of-safety based on ratio of medians
pos <- wilcox_test(I(log(asat)) ~ group, data = asat, alternative = "less", 
                   conf.int = TRUE, distribution = "exact")

### one-sided confidence set. Safety cannot be concluded since the effect of
### the compound exceeds 20% of the control median
exp(confint(pos)$conf.int)




cleanEx(); ..nameEx <- "expectation-methods"

### * expectation-methods

flush(stderr()); flush(stdout())

### Name: expectation-methods
### Title: Extract the Expectation and Covariance of Linear Statistics
### Aliases: expectation expectation-methods
###   expectation,IndependenceTest-method
###   expectation,ScalarIndependenceTestStatistic-method
###   expectation,MaxTypeIndependenceTestStatistic-method
###   expectation,QuadTypeIndependenceTestStatistic-method covariance
###   covariance-methods covariance,IndependenceTest-method
###   covariance,ScalarIndependenceTestStatistic-method
###   covariance,MaxTypeIndependenceTestStatistic-method
###   covariance,QuadTypeIndependenceTestStatistic-method
### Keywords: methods

### ** Examples


df <- data.frame(y = gl(3, 2), x = gl(3, 2)[sample(1:6)])

### Cochran-Mantel-Haenzel Test
ct <- cmh_test(y ~ x, data = df)
 
### the linear statistic, i.e, the contingency table
l <- statistic(ct, type = "linear")
l

### expectation
El <- expectation(ct)
El

### covariance
Vl <- covariance(ct)
Vl

### the standardized contingency table (hard way)
(l - El) / sqrt(variance(ct))

### easy way
statistic(ct, type = "standardized")




cleanEx(); ..nameEx <- "glioma"

### * glioma

flush(stderr()); flush(stdout())

### Name: glioma
### Title: Malignant Glioma Pilot Study
### Aliases: glioma
### Keywords: datasets

### ** Examples


data(glioma, package = "coin")

par(mfrow=c(1,2))

### Grade III glioma
g3 <- subset(glioma, histology == "Grade3")

### Plot Kaplan-Meier curves
plot(survfit(Surv(time, event) ~ group, data=g3), 
     main="Grade III Glioma", lty=c(2,1), 
     legend.text=c("Control", "Treated"),
     legend.bty=1, ylab="Probability", 
     xlab="Survival Time in Month")

### logrank test
surv_test(Surv(time, event) ~ group, data = g3, 
             distribution = "exact")

### Grade IV glioma
gbm <- subset(glioma, histology == "GBM")

### Plot Kaplan-Meier curves
plot(survfit(Surv(time, event) ~ group, data=gbm), 
     main="Grade IV Glioma", lty=c(2,1), 
     legend.text=c("Control", "Treated"),
     legend.bty=1, legend.pos=1, ylab="Probability", 
     xlab="Survival Time in Month")
   
### logrank test
surv_test(Surv(time, event) ~ group, data = gbm, 
             distribution = "exact")

### stratified logrank test
surv_test(Surv(time, event) ~ group | histology, data = glioma,
             distribution = approximate(B = 10000))




graphics::par(get("par.postscript", env = .CheckExEnv))
cleanEx(); ..nameEx <- "jobsatisfaction"

### * jobsatisfaction

flush(stderr()); flush(stdout())

### Name: jobsatisfaction
### Title: Income and Job Satisfaction
### Aliases: jobsatisfaction
### Keywords: datasets

### ** Examples


data(jobsatisfaction, package = "coin")

# Generalized Cochran-Mantel-Haenzel test
cmh_test(jobsatisfaction)




cleanEx(); ..nameEx <- "neuropathy"

### * neuropathy

flush(stderr()); flush(stdout())

### Name: neuropathy
### Title: Acute Painful Diabetic Neuropathy
### Aliases: neuropathy
### Keywords: datasets

### ** Examples


data(neuropathy, package = "coin")

### compare with Table 2 of Conover & Salsburg (1988)
oneway_test(pain ~ group, data = neuropathy, alternative = "less",
            distribution = "exact")

wilcox_test(pain ~ group, data = neuropathy, alternative = "less", 
            distribution = "exact")

oneway_test(pain ~ group, data = neuropathy, 
            distribution = approximate(B = 10000),
            alternative = "less", ytrafo = function(data) trafo(data,
            numeric_trafo = consal_trafo))




cleanEx(); ..nameEx <- "ocarcinoma"

### * ocarcinoma

flush(stderr()); flush(stdout())

### Name: ocarcinoma
### Title: Ovarian Carcinoma
### Aliases: ocarcinoma
### Keywords: datasets

### ** Examples


data(ocarcinoma, package = "coin")

### logrank test with exact two-sided p-value
lrt <- surv_test(Surv(time, event) ~ stadium, data = ocarcinoma,
                    distribution = "exact")

### the test statistic
statistic(lrt)

### p-value
pvalue(lrt)




cleanEx(); ..nameEx <- "pvalue-methods"

### * pvalue-methods

flush(stderr()); flush(stdout())

### Name: pvalue-methods
### Title: Extract P-Values
### Aliases: pvalue pvalue-methods pvalue,NullDistribution-method
###   pvalue,IndependenceTest-method pvalue,ScalarIndependenceTest-method
###   pvalue,MaxTypeIndependenceTest-method
###   pvalue,QuadTypeIndependenceTest-method
### Keywords: methods htest

### ** Examples


### artificial 2-sample problem
df <- data.frame(y = rnorm(20), x = gl(2, 10))
 
### Ansari-Bradley test
at <- ansari_test(y ~ x, data = df, distribution = "exact")

at

pvalue(at)




cleanEx(); ..nameEx <- "rotarod"

### * rotarod

flush(stderr()); flush(stdout())

### Name: rotarod
### Title: Rotating Rats Data
### Aliases: rotarod
### Keywords: datasets

### ** Examples


data(rotarod, package = "coin")

### Wilcoxon-Mann-Whitney Rank Sum Test

### one-sided exact (0.0186)
wilcox_test(time ~ group, data = rotarod, 
    alternative = "greater", distribution = "exact")
### two-sided exact (0.0373)
wilcox_test(time ~ group, data = rotarod, distribution = "exact")
### two-sided asymptotical (0.0147)
wilcox_test(time ~ group, data = rotarod)




cleanEx(); ..nameEx <- "sphase"

### * sphase

flush(stderr()); flush(stdout())

### Name: sphase
### Title: S-phase Fraction of Tumor Cells
### Aliases: sphase
### Keywords: datasets

### ** Examples


data(sphase, package = "coin")

maxstat_test(Surv(RFS, event) ~ SPF, data = sphase)




cleanEx(); ..nameEx <- "statistic-methods"

### * statistic-methods

flush(stderr()); flush(stdout())

### Name: statistic-methods
### Title: Extract Test Statistics, Linear Statistics and Standardized
###   Statistics
### Aliases: statistic statistic-methods statistic,IndependenceTest-method
###   statistic,ScalarIndependenceTestStatistic-method
###   statistic,MaxTypeIndependenceTestStatistic-method
###   statistic,QuadTypeIndependenceTestStatistic-method
### Keywords: methods

### ** Examples


df <- data.frame(y = gl(4, 5), x = gl(5, 4))

### Cochran-Mantel-Haenzel Test
ct <- cmh_test(y ~ x, data = df)

### chisq-type statistics
statistic(ct)

### the linear statistic, i.e, the contingency table
statistic(ct, type = "linear")

### the same
table(df$x, df$y)

### and the standardized contingency table for illustrating
### departures from the null-hypothesis of independence of x and y
statistic(ct, type = "standardized")



cleanEx(); ..nameEx <- "treepipit"

### * treepipit

flush(stderr()); flush(stdout())

### Name: treepipit
### Title: Tree Pipit (Anthus trivialis) Forest Data
### Aliases: treepipit
### Keywords: datasets

### ** Examples


data(treepipit, package = "coin")

maxstat_test(counts ~ age + coverstorey + coverregen + meanregen +
                      coniferous + deadtree + cbpiles + ivytree,
             data = treepipit)




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", env = .CheckExEnv),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
