\name{malformations}
\docType{data}
\alias{malformations}
\title{Maternal Drinking and Congenital Sex Organ Malformation}
\description{
  A subset of data from a study on the relationship between maternal alcohol
  consumption and congenital malformations.
}
\usage{malformations}
\format{
  A data frame with 32574 observations on 2 variables.
  \describe{
    \item{\code{consumption}}{
      alcohol consumption, an ordered factor with levels \code{"0"},
      \code{"<1"}, \code{"1-2"}, \code{"3-5"} and \code{">=6"}.
    }
    \item{\code{malformation}}{
      congenital sex organ malformation, a factor with levels \code{"Present"}
      and \code{"Absent"}.
    }
  }
}
\details{
  Data from a prospective study undertaken to determine whether moderate or
  light drinking during the first trimester of pregnancy increases the risk for
  congenital malformations (Mills and Graubard, 1987).  The subset given here
  concerns only sex organ malformation (Mills and Graubard, 1987, Tab. 4).

  Graubard and Korn (1987) used this data set to illustrate that different
  choices of scores for ordinal variables can lead to conflicting conclusions.
}
\source{
  Mills, J. L. and Graubard, B. I.  (1987).  Is moderate drinking during
  pregnancy associated with an increased risk for malformations?
  \emph{Pediatrics} \bold{80}(3), 309--314.
}
\references{
  Graubard, B. I. and Korn, E. L.  (1987).  Choice of column scores for testing
  independence in ordered \eqn{2 \times K}{2 x K} contingency tables.
  \emph{Biometrics} \bold{43}(2), 471--476.
}
\examples{
## Graubard and Korn (1987, Tab. 3)

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: midpoint scores (p < 0.05)
midpoints <- c(0, 0.5, 1.5, 4.0, 7.0)
chisq_test(malformation ~ consumption, data = malformations,
           distribution = approximate(B = 1000), alternative = "greater",
           scores = list(consumption = midpoints))

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: midrank scores (p > 0.05)
midranks <- c(8557.5, 24375.5, 32013.0, 32473.0, 32555.5)
chisq_test(malformation ~ consumption, data = malformations,
           distribution = approximate(B = 1000), alternative = "greater",
           scores = list(consumption = midranks))

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: equally spaced scores (p > 0.05)
chisq_test(malformation ~ consumption, data = malformations,
           distribution = approximate(B = 1000), alternative = "greater")
}
\keyword{datasets}
