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
  The data published by \bibcitet{coin::Hohnloser_1987} were 
  used by \bibcitet{coin::lausen_1992} to illustrate the use of
  maximally selected statistics.
}
\references{\bibshow{*}}
\examples{
## Asymptotic maximally selected logrank statistics
maxstat_test(Surv(time, event) ~ EF, data = hohnloser)
}
\keyword{datasets}
