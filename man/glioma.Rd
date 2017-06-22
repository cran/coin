\name{glioma}
\docType{data}
\alias{glioma}
\title{Malignant Glioma Pilot Study}
\description{
  A non-randomized pilot study on malignant glioma patients with pretargeted
  adjuvant radioimmunotherapy using yttrium-90-biotin.
}
\usage{glioma}
\format{
  A data frame with 37 observations on 7 variables.
  \describe{
    \item{\code{no.}}{
      patient number.
    }
    \item{\code{age}}{
      patient age (years).
    }
    \item{\code{sex}}{
      a factor with levels \code{"F"} (Female) and \code{"M"} (Male).
    }
    \item{\code{histology}}{
      a factor with levels \code{"GBM"} (grade IV) and \code{"Grade3"} (grade
      III).
    }
    \item{\code{group}}{
      a factor with levels \code{"Control"} and \code{"RIT"}.
    }
    \item{\code{event}}{
      status indicator for \code{time}: \code{FALSE} for right-censored
      observations and \code{TRUE} otherwise.
    }
    \item{\code{time}}{
      survival time (months).
    }
  }
}
\details{
  The primary endpoint of this small pilot study is survival.  Since the
  survival times are tied, the classical asymptotic logrank test may be
  inadequate in this setup.  Therefore, a permutation test using Monte Carlo
  resampling was computed in the original paper.  The data are taken from Tables
  1 and 2 of Grana \emph{et al.} (2002).
}
\source{
  Grana, C., Chinol, M., Robertson, C., Mazzetta, C., Bartolomei, M., De Cicco,
  C., Fiorenza, M., Gatti, M., Caliceti, P. and Paganelli, G.  (2002).
  Pretargeted adjuvant radioimmunotherapy with Yttrium-90-biotin in malignant
  glioma patients: A pilot study.  \emph{British Journal of Cancer}
  \bold{86}(2), 207--212.
}
\examples{
## Grade III glioma
g3 <- subset(glioma, histology == "Grade3")

## Plot Kaplan-Meier estimates
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:2, ncol = 2))
plot(survfit(Surv(time, event) ~ group, data = g3),
     main = "Grade III Glioma", lty = 2:1,
     ylab = "Probability", xlab = "Survival Time in Month",
     xlim = c(-2, 72))
legend("bottomleft", lty = 2:1, c("Control", "Treated"), bty = "n")

## Exact logrank test
logrank_test(Surv(time, event) ~ group, data = g3,
             distribution = "exact")


## Grade IV glioma
gbm <- subset(glioma, histology == "GBM")

## Plot Kaplan-Meier estimates
plot(survfit(Surv(time, event) ~ group, data = gbm),
     main = "Grade IV Glioma", lty = 2:1,
     ylab = "Probability", xlab = "Survival Time in Month",
     xlim = c(-2, 72))
legend("topright", lty = 2:1, c("Control", "Treated"), bty = "n")
par(op) # reset

## Exact logrank test
logrank_test(Surv(time, event) ~ group, data = gbm,
             distribution = "exact")


## Stratified approximative (Monte Carlo) logrank test
logrank_test(Surv(time, event) ~ group | histology, data = glioma,
             distribution = approximate(B = 10000))
}
\keyword{datasets}
