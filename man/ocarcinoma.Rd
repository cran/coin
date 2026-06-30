\name{ocarcinoma}
\docType{data}
\alias{ocarcinoma}
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
  Data from \bibcitet{coin::Fleming_1984} and \bibcitet{coin::Fleming_1980}
  were reanalysed by \bibcitet{coin::SchumacherSchulgen2002}.
}
\references{\bibshow{*}}
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
