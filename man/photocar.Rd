\name{photocar}
\alias{photocar}
\docType{data}
\title{ Multiple Dosing Photococarcinogenicity Experiment }
\description{
    Survival time, time to first tumor and total number of tumors for three
    groups of animals from a photococarcinogenicity study.
}
\usage{data("photocar")}
\format{
  A data frame with 108 observations on the following 6 variables.
  \describe{
    \item{\code{group}}{a factor with levels \code{A}, \code{B}, and \code{C}}
    \item{\code{ntumor}}{total number of tumors.}
    \item{\code{time}}{survival time.}
    \item{\code{event}}{censoring indicator (\code{TRUE} when the animal died).}
    \item{\code{dmin}}{time to first tumor.}
    \item{\code{tumor}}{censoring indicator for \code{dmin}, i.e.,
                        \code{TRUE} when at least one tumor was observed.}
  }
}
\details{

  The animals were exposed to different levels of ultraviolet radiation
  (UVR) exposure (group A:
  topical vehicle and 600 Robertson--Berger units of UVR, group B: no topical vehicle and 600
  Robertson--Berger units of UVR and group C: no topical vehicle and 1200 
  Robertson--Berger units of UVR). The
  data are taken from Tables 1-3 in Molefe et al. (2005).

  The main interest is testing the global null of no treatment effect with
  respect to survival time, time to first tumor and number of tumors (Molefe
  et al., 2005, analyse the detection time of tumors in addition, this data
  is not given here). In case the global null hypothesis can be rejected,
  the deviations from the partial hypotheses are of special interest.

}
\source{

    Daniel F. Molefe, James J. Chen, Paul C. Howard, Barbara J. Miller,
    Christopher P. Sambuco, P. Donald Forbes \& Ralph L. Kodell (2005).
    Tests for effects on tumor frequency and latency in multiple dosing
    photococarcinogenicity experiments. \emph{Journal of Statistical Planning
    and Inference} \bold{129}, 39--58.

    Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel \& Achim Zeileis (2006).
    A Lego system for conditional inference, \emph{The American Statistician},
    \bold{60}(3), 257--263.

}
\examples{

data("photocar")

layout(matrix(1:3, ncol = 3))
plot(survfit(Surv(time, event) ~ group, data = photocar), xmax = 50, lty =
1:3, main = "Survival Time")
legend("bottomleft", lty = 1:3, levels(photocar$group), bty = "n")
plot(survfit(Surv(dmin, tumor) ~ group, data = photocar), xmax = 50, lty =
1:3, main = "Time to First Tumor")
legend("bottomleft", lty = 1:3, levels(photocar$group), bty = "n")
boxplot(ntumor ~ group, data = photocar, main = "Number of Tumors")

### global test (all three responses)
fm <- Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group
it <- independence_test(fm, data = photocar, 
                        distribution = approximate(B = 10000))
pvalue(it)

### why was the global null hypothesis rejected?
statistic(it, "standardized")
pvalue(it, "single-step")

}
\keyword{datasets}
