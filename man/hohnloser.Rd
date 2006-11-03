\name{hohnloser}
\alias{hohnloser}
\title{
Left ventricular ejection fraction of patients with malignant ventricular
tachyarrhythmias.
}
\description{
A data frame with the left ventricular ejection fraction of patients with
malignant ventricular tachyarrhythmias including recurrence-free month and censoring.
}
\usage{
data("hohnloser")
}
\format{
\describe{
\item{\code{EF}}{left ventricular ejection in percent}
\item{\code{time}}{recurrence-free month}
\item{\code{event}}{censoring: 0 cencored, 1 not censored}
}
The data used here is published in Table 1 of Lausen and Schumacher (1992).
}
\source{

The data was first published by Hohnloser et al. (1987), 
the data used here is published in Table 1 of Lausen and Schumacher (1992).

}
\references{

    S. H. Hohnloser, E. A. Raeder, P. J. Podrid, T. B. Graboys \& B. and Lown,
    (1987). Predictors of antiarrhythmic drug efficacy in patients with
    malignant ventricular tachyarrhythmias. \emph{American Heart Journal}
    \bold{114}, 1--7.

    Berthold Lausen \& Martin Schumacher (1992). Maximally Selected Rank Statistics.
    \emph{Biometrics} \bold{48}, 73--85.

}

\examples{

  data("hohnloser", package = "coin")
  maxstat_test(Surv(time, event) ~ EF, data = hohnloser)

}
\keyword{datasets}
