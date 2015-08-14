\name{GTSG}
\docType{data}
\alias{GTSG}
\encoding{UTF-8}
\title{Gastrointestinal Tumor Study Group}
\description{
  A randomized clinical trial in gastric cancer.
}
\usage{GTSG}
\format{
  A data frame with 90 observations on 3 variables.
  \describe{
    \item{\code{time}}{
      survival time (days).
    }
    \item{\code{event}}{
      status indicator for \code{time}: \code{0} for censored observations and
      \code{1} otherwise.
    }
    \item{\code{group}}{
      a factor with levels \code{"Chemotherapy+Radiation"} and
      \code{"Chemotherapy"}.
    }
  }
}
\details{
  A clinical trial comparing chemotherapy alone versus a combination of
  chemotherapy and radiation therapy in the treatment of locally advanced,
  nonresectable gastric carcinoma.
}
\note{
  There is substantial separation between the estimated survival distributions
  at 8 to 10 months, but by month 26 the distributions intersect.
}
\source{
  Stablein, D. M., Carter, W. H., Jr. and Novak, J. W.  (1981).  Analysis of
  survival data with nonproportional hazard functions.  \emph{Controlled
    Clinical Trials} \bold{2}(2), 149--159.
}
\references{
  Moreau, T., Maccario, J., Lellouch, J. and Huber, C.  (1992).  Weighted log
  rank statistics for comparing two distributions.  \emph{Biometrika}
  \bold{79}(1), 195--198.

  Shen, W. and Le, C. T.  (2000).  Linear rank tests for censored survival data.
  \emph{Communications in Statistics -- Simulation and Computation}
  \bold{29}(1), 21--36.

  Tubert-Bitter, P., Kramar, A., \enc{Chalé}{Chale}, J. J. and Moureau, T.
  (1994).  Linear rank tests for comparing survival in two groups with crossing
  hazards.  \emph{Computational Statistics & Data Analysis} \bold{18}(5),
  547--559.
}
\examples{
## Plot Kaplan-Meier estimates
plot(survfit(Surv(time / (365.25 / 12), event) ~ group, data = GTSG),
     lty = 1:2, ylab = "\% Survival", xlab = "Survival Time in Months")
legend("topright", lty = 1:2,
       c("Chemotherapy+Radiation", "Chemotherapy"), bty = "n")

## Asymptotic logrank test
logrank_test(Surv(time, event) ~ group, data = GTSG)

## Asymptotic Prentice test
logrank_test(Surv(time, event) ~ group, data = GTSG, type = "Prentice")

## Asymptotic test against Weibull-type alternatives (Moreau et al., 1992)
moreau_weight <- function(time, n.risk, n.event)
    1 + log(-log(cumprod(n.risk / (n.risk + n.event))))

independence_test(Surv(time, event) ~ group, data = GTSG,
                  ytrafo = function(data)
                      trafo(data, surv_trafo = function(y)
                          logrank_trafo(y, weight = moreau_weight)))

## Asymptotic test against crossing-curve alternatives (Shen and Le, 2000)
shen_trafo <- function(x)
    ansari_trafo(logrank_trafo(x, type = "Prentice"))

independence_test(Surv(time, event) ~ group, data = GTSG,
                  ytrafo = function(data)
                      trafo(data, surv_trafo = shen_trafo))
}
\keyword{datasets}
