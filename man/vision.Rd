\name{vision}
\docType{data}
\alias{vision}
\encoding{UTF-8}
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
  1943--46, as given by Stuart (1953).

  This data set was used by Stuart (1955) to illustrate a test of marginal
  homogeneity.  Winell and \enc{Lindbäck}{Lindbaeck} (2018) also used the data,
  demonstrating a score-independent test for ordered categorical data.
}
\source{
  Stuart, A.  (1953).  The estimation and comparison of strengths of association
  in contingency tables.  \emph{Biometrika} \bold{40}(1/2), 105--110.
  \doi{10.2307/2333101}
}
\references{
  Stuart, A.  (1955).  A test for homogeneity of the marginal distributions in a
  two-way classification.  \emph{Biometrika} \bold{42}(3/4), 412--416.
  \doi{10.1093/biomet/42.3-4.412}

  Winell, H. and \enc{Lindbäck}{Lindbaeck}, J.  (2018).  A general
  score-independent test for order-restricted inference.  \emph{Statistics in
    Medicine} \bold{37}(21), 3078--3090.  \doi{10.1002/sim.7690}
}
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
