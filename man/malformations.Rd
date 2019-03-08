\name{malformations}
\docType{data}
\alias{malformations}
\encoding{UTF-8}
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

  This data set was used by Graubard and Korn (1987) to illustrate that
  different choices of scores for ordinal variables can lead to conflicting
  conclusions.  Zheng (2008) also used the data, demonstrating two different
  score-independent tests for ordered categorical data; see also Winell and
  \enc{Lindbäck}{Lindbaeck} (2018).
}
\source{
  Mills, J. L. and Graubard, B. I.  (1987).  Is moderate drinking during
  pregnancy associated with an increased risk for malformations?
  \emph{Pediatrics} \bold{80}(3), 309--314.
}
\references{
  Graubard, B. I. and Korn, E. L.  (1987).  Choice of column scores for testing
  independence in ordered \eqn{2 \times K}{2 x K} contingency tables.
  \emph{Biometrics} \bold{43}(2), 471--476.  \doi{10.2307/2531828}

  Winell, H. and \enc{Lindbäck}{Lindbaeck}, J.  (2018).  A general
  score-independent test for order-restricted inference.  \emph{Statistics in
    Medicine} \bold{37}(21), 3078--3090.  \doi{10.1002/sim.7690}

  Zheng, G.  (2008).  Analysis of ordered categorical data: Two
  score-independent approaches.  \emph{Biometrics} \bold{64}(4), 1276–-1279.
  \doi{10.1111/j.1541-0420.2008.00992.x}
}
\examples{
## Graubard and Korn (1987, Tab. 3)

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: midpoint scores (p < 0.05)
midpoints <- c(0, 0.5, 1.5, 4.0, 7.0)
chisq_test(malformation ~ consumption, data = malformations,
           distribution = approximate(nresample = 1000),
           alternative = "greater",
           scores = list(consumption = midpoints))

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: midrank scores (p > 0.05)
midranks <- c(8557.5, 24375.5, 32013.0, 32473.0, 32555.5)
chisq_test(malformation ~ consumption, data = malformations,
           distribution = approximate(nresample = 1000),
           alternative = "greater",
           scores = list(consumption = midranks))

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: equally spaced scores (p > 0.05)
chisq_test(malformation ~ consumption, data = malformations,
           distribution = approximate(nresample = 1000),
           alternative = "greater")

\dontrun{
## One-sided approximative (Monte Carlo) score-independent test
## Winell and Lindbaeck (2018)
(it <- independence_test(malformation ~ consumption, data = malformations,
                         distribution = approximate(nresample = 1000,
                                                    parallel = "snow",
                                                    ncpus = 8),
                         alternative = "greater",
                         xtrafo = function(data)
                             trafo(data, ordered_trafo = zheng_trafo)))

## Extract the "best" set of scores
ss <- statistic(it, type = "standardized")
idx <- which(ss == max(ss), arr.ind = TRUE)
ss[idx[1], idx[2], drop = FALSE]}
}
\keyword{datasets}
