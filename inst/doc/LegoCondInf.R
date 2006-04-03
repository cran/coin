###################################################
### chunk number 1: setup
###################################################
options(width = 65, prompt = "R> ")
require("coin")
set.seed(290875)
### get rid of the NAMESPACE
load(file.path(.find.package("coin"), "R", "all.rda"))
anonymous <- FALSE


###################################################
### chunk number 2: authors
###################################################
if(!anonymous) {
    cat("\\author{Torsten Hothorn$^1$, Kurt Hornik$^2$, \\\\ 
            Mark A. van de Wiel$^3$ and Achim Zeileis$^2$}\n")
} else  {
    cat("\\author{Anonymous}\n")
}


###################################################
### chunk number 3: affil
###################################################
if(!anonymous)
    cat("\\noindent$^1$ Institut f\\\"ur Medizininformatik, Biometrie und Epidemiologie\\\\
           Friedrich-Alexander-Universit\\\"at Erlangen-N\\\"urnberg\\\\
           Waldstra{\\ss}e 6, D-91054 Erlangen, Germany \\\\
           \\texttt{Torsten.Hothorn@R-project.org}
         \\newline

         \\noindent$^2$ Department f\\\"ur Statistik und Mathematik,
            Wirtschaftsuniversit\\\"at Wien \\\\
            Augasse 2-6, A-1090 Wien, Austria \\\\
            \\texttt{Kurt.Hornik@R-project.org} \\\\
            \texttt{Achim.Zeileis@R-project.org}
         \\newline

         \\noindent$^3$ Department of Mathematics and Computer Science \\\\
            Eindhoven University of Technology \\\\
            HG 9.25, P.O. Box 513 \\\\
            5600 MB Eindhoven, The Netherlands \\\\
            \\texttt{m.a.v.d.wiel@TUE.nl}
         \\newline\n")


###################################################
### chunk number 4: coincite
###################################################
if (anonymous) {
    cat(" \\citep{PKG:coina} ")
} else {
    cat(" \\citep{PKG:coin} ")
}


###################################################
### chunk number 5: alpha-data-figure
###################################################
n <- table(alpha$alength)
par(cex.lab = 1.3, cex.axis = 1.3)
boxplot(elevel ~ alength, data = alpha, ylab = "Expression Level",
        xlab = "NACP-REP1 Allele Length", varwidth = TRUE)
axis(3, at = 1:3, labels = paste("n = ", n))
rankif <- function(data) trafo(data, numeric_trafo = rank)


###################################################
### chunk number 6: alpha-kruskal
###################################################
independence_test(elevel ~ alength, data = alpha, ytrafo = function(data)
                  trafo(data, numeric_trafo = rank), teststat = "quadtype")


###################################################
### chunk number 7: alpha-kruskal
###################################################
kruskal.test(elevel ~ alength, data = alpha)


###################################################
### chunk number 8: alpha-kruskal-ordered
###################################################
independence_test(elevel ~ alength, data = alpha, ytrafo = function(data)
                  trafo(data, numeric_trafo = rank), 
                  scores = list(alength = c(2, 7, 11)))


###################################################
### chunk number 9: alzheimer-demographics
###################################################
total <- nrow(alzheimer)
stopifnot(total == 538)
male <- sum(alzheimer$gender == "Male")
stopifnot(male == 200)
female <- sum(alzheimer$gender == "Female")
stopifnot(female == 338)
disease <- table(alzheimer$disease)
smoked <- sum(alzheimer$smoking != "None")
atab <- xtabs(~ smoking + + disease + gender, data = alzheimer)
### there is a discrepancy between Table 1 (32% smokers of 117 women 
### suffering from other diagnoses) and Table 4 (63% non-smokers).
### We used the data as given in Table 4.


###################################################
### chunk number 10: alzheimer-tab
###################################################
x <- t(atab[,,"Female"])
lines <- paste(paste(dimnames(x)$disease, " & "), 
               paste(apply(x, 1, function(l) paste(l, collapse = " & ")), "\\\\"))
for (i in 1:length(lines)) cat(lines[i], "\n")


###################################################
### chunk number 11: alzheimer-tab
###################################################
x <- t(atab[,,"Male"])
lines <- paste(paste(dimnames(x)$disease, " & "), 
               paste(apply(x, 1, function(l) paste(l, collapse = " & ")), "\\\\"))
for (i in 1:length(lines)) cat(lines[i], "\n")


###################################################
### chunk number 12: alzheimer-plot
###################################################
layout(matrix(1:2, ncol = 2))
spineplot(disease ~ smoking, data = alzheimer, subset = gender == "Male",
main = "Male", xlab = "Smoking", ylab = "Disease", tol = 1)
spineplot(disease ~ smoking, data = alzheimer, subset = gender == "Female",
main = "Female", xlab = "Smoking", ylab = "Disease", tol = 1)


###################################################
### chunk number 13: alzheimer-mantelhaen
###################################################
it_alz <- independence_test(disease ~ smoking | gender, data = alzheimer, 
                            teststat = "quadtype")
it_alz


###################################################
### chunk number 14: alzheimer-statistic
###################################################
statistic(it_alz, type = "linear")


###################################################
### chunk number 15: alzheimer-men
###################################################
females <- alzheimer$gender == "Female"
pvalue(independence_test(disease ~ smoking, data = alzheimer,
       subset = females, teststat = "quadtype"))
pvalue(independence_test(disease ~ smoking, data = alzheimer,
       subset = !females, teststat = "quadtype"))


###################################################
### chunk number 16: alzheimer-max
###################################################
it_alzmax <- independence_test(disease ~ smoking, data = alzheimer,
       subset = !females, teststat = "maxtype")
it_alzmax


###################################################
### chunk number 17: alzheimer-maxstat
###################################################
statistic(it_alzmax, "standardized")


###################################################
### chunk number 18: alzheimer-qperm
###################################################
qperm(it_alzmax, 0.95)


###################################################
### chunk number 19: alzheimer-MTP
###################################################
pvalue(it_alzmax, method = "single-step")


###################################################
### chunk number 20: photocar-plot
###################################################
par(cex.lab = 1.3, cex.axis = 1.3)
layout(matrix(1:3, ncol = 3))
plot(survfit(Surv(time, event) ~ group, data = photocar), xmax = 50,
     xlab = "Survival Time (in weeks)", ylab = "Probability",
     lty = 1:3)
     legend("bottomleft", lty = 1:3, levels(photocar$group), bty = "n")
plot(survfit(Surv(dmin, tumor) ~ group, data = photocar), xmax = 50,
     xlab = "Time to First Tumor (in weeks)", ylab = "Probability",
     lty = 1:3)
     legend("bottomleft", lty = 1:3, levels(photocar$group), bty = "n")
boxplot(ntumor ~ group, data = photocar, 
        ylab = "Number of Tumors", xlab = "Treatment Group")


###################################################
### chunk number 21: photocar-global
###################################################
it_ph <- independence_test(Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group,
                           data = photocar)
it_ph 


###################################################
### chunk number 22: photocar-linear
###################################################
statistic(it_ph, type = "linear")


###################################################
### chunk number 23: photocar-stand
###################################################
statistic(it_ph, type = "standardized")


###################################################
### chunk number 24: photocar-stand eval=FALSE
###################################################
## pvalue(it_ph, method = "single-step")


###################################################
### chunk number 25: photocar-stand
###################################################
round(pvalue(it_ph, method = "single-step"), 5)


###################################################
### chunk number 26: photocar-MC
###################################################
it <- independence_test(Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group,
                  data = photocar, distribution = approximate(50000))
pvalue(it, method = "single-step")


###################################################
### chunk number 27: photocar-MC2
###################################################
pvalue(it, method = "step-down")


###################################################
### chunk number 28: mercuryfish-plot
###################################################
par(cex.lab = 1.3, cex.axis = 1.3)
layout(matrix(1:3, ncol = 3))
boxplot(I(log(mercury)) ~ group, data = mercuryfish, 
        ylab = "Mercury Blood Level (log scale)")
boxplot(abnormal ~ group, data = mercuryfish, 
        ylab = "Abnormal Cells (in %)")
boxplot(ccells ~ group, data = mercuryfish, 
        ylab = "Chromosome Aberrations (in %)")


###################################################
### chunk number 29: mercurysfish-score
###################################################
coherence <- function(data) {
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


###################################################
### chunk number 32: mercuryfish-ppermplot
###################################################
par(cex.lab = 1.3, cex.axis = 1.1)
ite <- poset
ita <- independence_test(mercury + abnormal + ccells ~ group, data =     
                           mercuryfish, ytrafo = coherence)
site <- support(ite)
layout(matrix(1:2, ncol = 2))
site <- site[site <= qperm(ite, 0.1) & site > -3]
pite <- sapply(site, function(x) pperm(ite, x))
pita <- sapply(site, function(x) pperm(ita, x))

plot(site, pite, type = "S", ylab = "Probability", xlab = "Standardized Statistic")
lines(site, pita, lty = 3)
legend("topleft", lty = c(1,3), legend = c("Conditional Distribution",
"Approximation"), bty = "n")

site <- support(ite)
site <- site[site >= qperm(ite, 0.9) & site < 3]
pite <- sapply(site, function(x) pperm(ite, x))
pita <- sapply(site, function(x) pperm(ita, x))

plot(site, pite, type = "S", ylab = "Probability", xlab = "Standardized Statistic")
lines(site, pita, lty = 3)


