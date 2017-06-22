\name{hohnloser}
\docType{data}
\alias{hohnloser}
\title{Left Ventricular Ejection Fraction}
\description{
  Left ventricular ejection fraction in patients with malignant ventricular
  tachyarrhythmias including recurrence-free month and censoring.
}
\usage{hohnloser}
\format{
  A data frame with 94 observations on 3 variables.
  \describe{
    \item{\code{EF}}{
      ejection fraction (\%).
    }
    \item{\code{time}}{
      recurrence-free month.
    }
    \item{\code{event}}{
      status indicator for \code{time}: \code{0} for right-censored observations
      and \code{1} otherwise.
    }
  }
}
\details{
  The data was used by Lausen and Schumacher (1992) to illustrate the use of
  maximally selected statistics.
}
\source{
  Hohnloser, S. H., Raeder, E. A., Podrid, P. J., Graboys, T. B. and Lown, B.
  (1987).  Predictors of antiarrhythmic drug efficacy in patients with malignant
  ventricular tachyarrhythmias.  \emph{American Heart Journal} \bold{114}(1 Pt
  1), 1--7.
}
\references{
  Lausen, B. and Schumacher, M.  (1992).  Maximally selected rank statistics.
  \emph{Biometrics} \bold{48}(1), 73--85.
}
\examples{
## Asymptotic maximally selected logrank statistics
maxstat_test(Surv(time, event) ~ EF, data = hohnloser)
}
\keyword{datasets}
