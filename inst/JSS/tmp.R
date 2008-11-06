###################################################
### chunk number 1: coin-setup
###################################################
options(prompt = "R> ", continue = "+  ")
library("coin")
library("e1071")
set.seed(290875)
### get rid of the NAMESPACE
#attach(asNamespace("coin")) ## FIXME: This must not be necessary!


###################################################
### chunk number 2: Ex
###################################################
library("coin")
data("rotarod", package = "coin")
independence_test(time ~ group, data = rotarod,
  ytrafo = rank, distribution = exact())


###################################################
### chunk number 3: IndependenceProblem
###################################################


###################################################
### chunk number 4: Ex-IndependenceProblem
###################################################
ip <- new("IndependenceProblem",
  y = rotarod["time"], x = rotarod["group"])


###################################################
### chunk number 5: IndependenceTestProblem
###################################################


###################################################
### chunk number 6: Ex-IndependenceTestProblem
###################################################
itp <- new("IndependenceTestProblem", ip, ytrafo = rank)


###################################################
### chunk number 7: IndependenceLinearStatistic
###################################################


###################################################
### chunk number 8: IndependenceTestStatistic
###################################################


###################################################
### chunk number 9: Ex-IndependenceTestStatistic
###################################################
its <- new("IndependenceTestStatistic", itp)


###################################################
### chunk number 10: Ex-IndependenceTestStatistic-statistic
###################################################
statistic(its, "linear")


###################################################
### chunk number 11: Ex-IndependenceTestStatistic-statistic
###################################################
expectation(its)
variance(its)


###################################################
### chunk number 12: ScalarIndependenceTestStatistic
###################################################


###################################################
### chunk number 13: Ex-ScalarIndependenceTestStatistic
###################################################
sits <- new("ScalarIndependenceTestStatistic", its,
  alternative = "two.sided")
statistic(sits, "standardized")


###################################################
### chunk number 14: MaxTypeIndependenceTestStatistic
###################################################


###################################################
### chunk number 15: QuadTypeIndependenceTestStatistic
###################################################


###################################################
### chunk number 16: PValue
###################################################


###################################################
### chunk number 17: NullDistribution
###################################################


###################################################
### chunk number 18: Ex-NullDistribution-pvalue
###################################################
end <- ExactNullDistribution(sits)
pvalue(end, statistic(sits))
qperm(end, 0.95)


###################################################
### chunk number 19: IndependenceTest
###################################################


###################################################
### chunk number 20: IndependenceTest
###################################################
new("IndependenceTest", statistic = sits, distribution = end)


###################################################
### chunk number 21: Ex-distribution
###################################################
set.seed(2908)
correxample <- data.frame(x = rnorm(7), y = rnorm(7))
sexact <- function(object) {
  x <- object@xtrans  
  y <- object@ytrans  
  perms <- permutations(nrow(x))
  pstats <- apply(perms, 1, function(p) sum(x[p,] * y))
  pstats <- (pstats - expectation(object)) / sqrt(variance(object)) 
  p <- function(q) 1 - mean(pstats > q)
  new("PValue", p = p, pvalue = p)
}


###################################################
### chunk number 22: Ex-distribution
###################################################
independence_test(y ~ x, data = correxample, alternative = "less",
  distribution = sexact)


###################################################
### chunk number 23: 
###################################################
mood_score <- function(y) (rank(y) - (nrow(y) + 1) / 2)^2


###################################################
### chunk number 24: 
###################################################
ip <- new("IndependenceProblem",
  y = rotarod["time"], x = rotarod["group"])
itp <- new("IndependenceTestProblem", ip, 
  ytrafo = mood_score)
its <- new("IndependenceTestStatistic", itp)
sits <- new("ScalarIndependenceTestStatistic", its, 
  alternative = "two.sided")
new("ScalarIndependenceTest", statistic = sits,
  distribution = ExactNullDistribution(sits, algorithm = "split-up"))


###################################################
### chunk number 25: 
###################################################
independence_test(time ~ group, data = rotarod, ytrafo = mood_score,
  distribution = exact(algorithm = "split-up"))


###################################################
### chunk number 26: js
###################################################
data("jobsatisfaction", package = "coin")
js <- jobsatisfaction
dimnames(js)[[2]] <- c("VeryDiss", "ModDiss", "ModSat", "VerySat")
ftable(Job.Satisfaction ~ Gender + Income, data = js)


###################################################
### chunk number 27: js-plot
###################################################
library("vcd")
cotabplot(js,
  split_vertical = TRUE, spacing = spacing_highlighting,
  gp = gpar(fill = rev(gray.colors(4))),
  labeling_args = list(rot_labels = 0, varnames = FALSE,
    just_labels = c("center", "right")),  
  panel_args = list(margins = c(3, 1, 2, 3.5)))


###################################################
### chunk number 28: jobsatisfaction-it
###################################################
it <- independence_test(js, teststat = "quad",
  distribution = asymptotic())
it


###################################################
### chunk number 29: jobsatisfaction-T
###################################################
statistic(it, "linear")


###################################################
### chunk number 30: jobsatisfaction-margin
###################################################
margin.table(js, 1:2)


###################################################
### chunk number 31: jobsatisfaction-stat
###################################################
statistic(it, "standardized")


###################################################
### chunk number 32: jobsatisfaction-max
###################################################
independence_test(js, teststat = "max")


###################################################
### chunk number 33: jobsatisfaction-minp
###################################################
pvalue(independence_test(js, teststat = "max"), method = "single-step")


###################################################
### chunk number 34: jobsatisfaction-ordinal
###################################################
it <- independence_test(js, distribution = approximate(B = 10000),
  scores = list(Job.Satisfaction = c(1, 3, 4, 5),
    Income = c(3, 10, 20, 35)))
pvalue(it)


###################################################
### chunk number 35: coin-doxygen eval=FALSE
###################################################
## browseURL(system.file("documentation", "html", "index.html",
##   package = "coin"))


