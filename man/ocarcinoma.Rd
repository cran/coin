\name{ocarcinoma}
\docType{data}
\alias{ocarcinoma}
\encoding{UTF-8}
\title{Ovarian Carcinoma}
\description{
  Survival times of 35 women suffering from ovarian carcinoma at stadium II and
  IIA.
}
\usage{ocarcinoma}
\format{
  A data frame with 35 observations on 3 variables.
  \describe{
    \item{\code{time}}{
      time (days).
    }
    \item{\code{stadium}}{
      a factor with levels \code{"II"} and \code{"IIA"}.
    }
    \item{\code{event}}{
      status indicator for \code{time}: \code{FALSE} for right-censored
      observations and \code{TRUE} otherwise.
    }
  }
}
\details{
  Data from Fleming \emph{et al.} (1980) and Fleming, Green and Harrington
  (1984).  Reanalysed in Schumacher and Schulgen (2002).
}
\source{
  Fleming, T. R., Green, S. J. and Harrington, D. P.  (1984).  Considerations
  for monitoring and evaluating treatment effects in clinical trials.
  \emph{Controlled Clinical Trials} \bold{5}(1), 55--66.
  \doi{10.1016/0197-2456(84)90150-8}

  Fleming, T. R., O'Fallon, J. R., O'Brien, P. C. and Harrington, D. P.  (1980).
  Modified Kolmogorov-Smirnov test procedures with application to arbitrarily
  right-censored data.  \emph{Biometrics} \bold{36}(4), 607--625.
  \doi{10.2307/2556114}
}
\references{
  Schumacher, M. and Schulgen, G.  (2002).  \emph{Methodik Klinischer Studien:
    Methodische Grundlagen der Planung, \enc{Durchführung}{Durchfuehrung} und
    Auswertung}.  Heidelberg: Springer.
}
\examples{
## Exact logrank test
lt <- logrank_test(Surv(time, event) ~ stadium, data = ocarcinoma,
                   distribution = "exact")

## Test statistic
statistic(lt)

## P-value
pvalue(lt)
}
\keyword{datasets}
