\name{mercuryfish}
\docType{data}
\alias{mercuryfish}
\title{Chromosomal Effects of Mercury-Contaminated Fish Consumption}
\description{
  The mercury level in blood, the proportion of cells with abnormalities, and
  the proportion of cells with chromosome aberrations in consumers of
  mercury-contaminated fish and a control group.
}
\usage{mercuryfish}
\format{
  A data frame with 39 observations on 4 variables.
  \describe{
    \item{\code{group}}{
      a factor with levels \code{"control"} and \code{"exposed"}.
    }
    \item{\code{mercury}}{
      mercury level in blood.
    }
    \item{\code{abnormal}}{
      the proportion of cells with structural abnormalities.
    }
    \item{\code{ccells}}{
      the proportion of \eqn{C_u} cells, i.e., cells with asymmetrical or
      incomplete-symmetrical chromosome aberrations.
    }
  }
}
\details{
  Control subjects (\code{"control"}) and subjects who ate contaminated fish for
  more than three years (\code{"exposed"}) are under study.

  Rosenbaum (1994) proposed a coherence criterion defining a partial ordering,
  i.e., an observation is smaller than another when all responses are smaller,
  and a score reflecting the \dQuote{ranking} is attached to each observation.
  The corresponding partially ordered set (POSET) test can be used to test if
  the distribution of the scores differ between the groups.  Alternatively, a
  multivariate test can be applied.
}
\source{
  Skerfving, S., Hansson, K., Mangs, C., Lindsten, J. and Ryman, N.  (1974).
  Methylmercury-induced chromosome damage in men.  \emph{Environmental Research}
  \bold{7}(1), 83--98.  \doi{10.1016/0013-9351(74)90078-4}
}
\references{
  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2006).  A Lego
  system for conditional inference.  \emph{The American Statistician}
  \bold{60}(3), 257--263.  \doi{10.1198/000313006X118430}

  Rosenbaum, P. R. (1994).  Coherence in observational studies.
  \emph{Biometrics} \bold{50}(2), 368--374.  \doi{10.2307/2533380}
}
\examples{
## Coherence criterion
coherence <- function(data) {
    x <- as.matrix(data)
    matrix(apply(x, 1, function(y)
        sum(colSums(t(x) < y) == ncol(x)) -
            sum(colSums(t(x) > y) == ncol(x))), ncol = 1)
}

## Asymptotic POSET test
poset <- independence_test(mercury + abnormal + ccells ~ group,
                           data = mercuryfish, ytrafo = coherence)

## Linear statistic (T in the notation of Rosenbaum, 1994)
statistic(poset, type = "linear")

## Expectation
expectation(poset)

## Variance
## Note: typo in Rosenbaum (1994, p. 371, Sec. 2, last paragraph)
variance(poset)

## Standardized statistic
statistic(poset)

## P-value
pvalue(poset)

## Exact POSET test
independence_test(mercury + abnormal + ccells ~ group,
                  data = mercuryfish, ytrafo = coherence,
                  distribution = "exact")

## Asymptotic multivariate test
mvtest <- independence_test(mercury + abnormal + ccells ~ group,
                            data = mercuryfish)

## Global p-value
pvalue(mvtest)

## Single-step adjusted p-values
pvalue(mvtest, method = "single-step")

## Step-down adjusted p-values
pvalue(mvtest, method = "step-down")
}
\keyword{datasets}
