###################################################
### chunk number 1: setup
###################################################
options(width = 65, prompt = "R> ")
require("coin")
set.seed(290875)
### get rid of the NAMESPACE
#load(file.path(.find.package("coin"), "R", "all.rda"))
#anonymous <- FALSE



###################################################
### chunk number 29: mercurysfish-score
###################################################
coherence <- function(data) {
    print(data)
    x <- t(as.matrix(data))
    f <- function(y) 
        sum(colSums(x < y) == nrow(x)) - sum(colSums(x > y) == nrow(x))
    matrix(apply(x, 2, f), ncol = 1)
}


###################################################
### chunk number 30: mercuryfish-poset
###################################################
poset <- independence_test(mercury + abnormal + ccells ~ group, 
    data = mercuryfish, ytrafo = coherence, distribution = exact())


###################################################
### chunk number 31: mercuryfish-pvalue
###################################################
pvalue(poset)

