\name{vision}
\docType{data}
\alias{vision}
\title{Unaided Distance Vision}
\description{
  Assessment of unaided distance vision of women in Britain.
}
\usage{vision}
\format{
  A contingency table with 7477 observations on 2 variables.
  \describe{
    \item{\code{Right.Eye}}{
      a factor with levels \code{"Highest Grade"}, \code{"Second Grade"},
      \code{"Third Grade"} and \code{"Lowest Grade"}.
    }
    \item{\code{Left.Eye}}{
      a factor with levels \code{"Highest Grade"},  \code{"Second Grade"},
      \code{"Third Grade"} and \code{"Lowest Grade"}.
    }
  }
}
\details{
  Paired ordered categorical data from case-records of eye-testing of 7477 women
  aged 30--39 years employed by Royal Ordnance Factories in Britain during
  1943--46, as given by \bibcitet{coin::stuart_1955}.

  This data set was used by  \bibcitet{coin::stuart_1955} to illustrate a test of marginal
  homogeneity.  \bibcitet{coin::winell_2018} also used the data,
  demonstrating a score-independent test for ordered categorical data.
}
\references{\bibshow{*}}
\examples{
## Asymptotic Stuart test (Q = 11.96)
diag(vision) <- 0 # speed-up
mh_test(vision)

## Asymptotic score-independent test
## Winell and Lindbaeck (2018)
(st <- symmetry_test(vision,
                     ytrafo = function(data)
                         trafo(data, factor_trafo = function(y)
                             zheng_trafo(as.ordered(y)))))
ss <- statistic(st, type = "standardized")
idx <- which(abs(ss) == max(abs(ss)), arr.ind = TRUE)
ss[idx[1], idx[2], drop = FALSE]
}
\keyword{datasets}
