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
  congenital malformations \bibcitep{coin::Mills+Graubard:1987}.  The subset given here
  concerns only sex organ malformation
  \bibcitep{|coin::Mills+Graubard:1987|Tab. 4}.

  This data set was used by \bibcitet{coin::graubard_1987} to illustrate that
  different choices of scores for ordinal variables can lead to conflicting
  conclusions.  \bibcitet{coin::zheng:2008} also used the data, demonstrating two different
  score-independent tests for ordered categorical data; see also
  \bibcitet{coin::winell_2018}.
}
\references{\bibshow{*}}
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
