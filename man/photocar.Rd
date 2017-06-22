\name{photocar}
\docType{data}
\alias{photocar}
\title{Multiple Dosing Photococarcinogenicity Experiment}
\description{
  Survival time, time to first tumor, and total number of tumors in three groups
  of animals in a photococarcinogenicity study.
}
\usage{photocar}
\format{
  A data frame with 108 observations on 6 variables.
  \describe{
    \item{\code{group}}{
      a factor with levels \code{"A"}, \code{"B"}, and \code{"C"}.
    }
    \item{\code{ntumor}}{
      total number of tumors.
    }
    \item{\code{time}}{
      survival time.
    }
    \item{\code{event}}{
      status indicator for \code{time}: \code{FALSE} for right-censored
      observations and \code{TRUE} otherwise.
    }
    \item{\code{dmin}}{
      time to first tumor.
    }
    \item{\code{tumor}}{
      status indicator for \code{dmin}: \code{FALSE} when no tumor was observed
      and \code{TRUE} otherwise.
    }
  }
}
\details{
  The animals were exposed to different levels of ultraviolet radiation (UVR)
  exposure (group A: topical vehicle and 600 Robertson--Berger units of UVR,
  group B: no topical vehicle and 600 Robertson--Berger units of UVR and group
  C: no topical vehicle and 1200 Robertson--Berger units of UVR).  The data are
  taken from Tables 1 to 3 in Molefe \emph{et al.} (2005).

  The main interest is testing the global null hypothesis of no treatment effect
  with respect to survival time, time to first tumor and number of tumors.
  (Molefe \emph{et al.}, 2005, also analysed the detection time of tumors, but
  that data is not given here.)  In case the global null hypothesis can be
  rejected, the deviations from the partial null hypotheses are of special
  interest.
}
\source{
  Molefe, D. F., Chen, J. J., Howard, P. C., Miller, B. J., Sambuco, C. P.,
  Forbes, P. D. and Kodell, R. L.  (2005).  Tests for effects on tumor frequency
  and latency in multiple dosing photococarcinogenicity experiments.
  \emph{Journal of Statistical Planning and Inference} \bold{129}(1--2), 39--58.
}
\references{
  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2006).  A Lego
  system for conditional inference.  \emph{The American Statistician}
  \bold{60}(3), 257--263.
}
\examples{
## Plotting data
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:3, ncol = 3))
with(photocar, {
    plot(survfit(Surv(time, event) ~ group),
         lty =  1:3, xmax = 50, main = "Survival Time")
    legend("bottomleft", lty = 1:3, levels(group), bty = "n")
    plot(survfit(Surv(dmin, tumor) ~ group),
         lty = 1:3, xmax = 50, main = "Time to First Tumor")
    legend("bottomleft", lty = 1:3, levels(group), bty = "n")
    boxplot(ntumor ~ group, main = "Number of Tumors")
})
par(op) # reset

## Approximative multivariate (all three responses) test
it <- independence_test(Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group,
                        data = photocar,
                        distribution = approximate(B = 10000))

## Global p-value
pvalue(it)

## Why was the global null hypothesis rejected?
statistic(it, type = "standardized")
pvalue(it, method = "single-step")
}
\keyword{datasets}
