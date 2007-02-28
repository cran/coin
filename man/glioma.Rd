\name{glioma}
\alias{glioma}
\docType{data}
\title{ Malignant Glioma Pilot Study}
\description{
  A non-randomized pilot study on malignant glioma patients with 
  pretargeted adjuvant radioimmunotherapy using Yttrium-90-biotin.
}
\usage{data("glioma")}
\format{
  A data frame with 37 observations on the following 7 variables.
  \describe{
    \item{no.}{patient number.}
    \item{age}{patients ages in years.}
    \item{sex}{a factor with levels \code{F}(emale) and \code{M}(ale). }
    \item{histology}{a factor with levels \code{GBM} (grade IV) and 
                     \code{Grade3} (grade III)}
    \item{time}{survival times in month.}
    \item{event}{censoring indicator: \code{FALSE} censored and \code{TRUE} dead.}
    \item{group}{a factor with levels \code{Control} and \code{RIT}.}
  }
}
\details{
  The primary endpoint of this small pilot study is survival. 
  Survival times are tied, the usual asymptotic log-rank test may be
  inadequate in this setup. Therefore, a permutation test (via Monte-Carlo sampling) 
  was conducted in the original paper. The data are taken from Tables 1 and 2 of 
  Grana et al. (2002).
}
\source{

    C. Grana, M. Chinol, C. Robertson, C. Mazzetta, M. Bartolomei, C. De
    Cicco, M. Fiorenza, M. Gatti, P. Caliceti \& G. Paganelli (2002),
    Pretargeted adjuvant radioimmunotherapy with Yttrium-90-biotin in malignant
    glioma patients: A pilot study. \emph{British Journal of Cancer}
    \bold{86}(2), 207--212.
}
\examples{

  layout(matrix(1:2, ncol = 2))

  ### Grade III glioma
  g3 <- subset(glioma, histology == "Grade3")

  ### Plot Kaplan-Meier curves
  plot(survfit(Surv(time, event) ~ group, data = g3), 
       main = "Grade III Glioma", lty = c(2,1), 
       legend.text = c("Control", "Treated"),
       legend.bty = 1, ylab = "Probability", 
       xlab = "Survival Time in Month")

  ### logrank test
  surv_test(Surv(time, event) ~ group, data = g3, 
               distribution = "exact")

  ### Grade IV glioma
  gbm <- subset(glioma, histology == "GBM")

  ### Plot Kaplan-Meier curves
  plot(survfit(Surv(time, event) ~ group, data = gbm), 
       main = "Grade IV Glioma", lty = c(2,1), 
       legend.text = c("Control", "Treated"),
       legend.bty = 1, legend.pos = 1, ylab = "Probability", 
       xlab = "Survival Time in Month")
   
  ### logrank test
  surv_test(Surv(time, event) ~ group, data = gbm, 
            distribution = "exact")

  ### stratified logrank test
  surv_test(Surv(time, event) ~ group | histology, data = glioma,
            distribution = approximate(B = 10000))

}
\keyword{datasets}
