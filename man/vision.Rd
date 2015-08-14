\name{vision}
\alias{vision}
\docType{data}
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

  Stuart (1955) used this data set to illustrate a test of marginal homogeneity.
}
\source{
  Stuart, A.  (1953).  The estimation and comparison of strengths of association
  in contingency tables.  \emph{Biometrika} \bold{40}(1/2), 105--110.
}
\references{
  Stuart, A.  (1955).  A test for homogeneity of the marginal distribution in a
  two-way classification.  \emph{Biometrika} \bold{42}(3/4), 412--416.
}
\examples{
## Asymptotic Stuart(-Maxwell) test (Q = 11.96)
diag(vision) <- 0 # speed-up
mh_test(vision)
}
\keyword{datasets}
