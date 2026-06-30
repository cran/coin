\name{alpha}
\docType{data}
\alias{alpha}
\title{Genetic Components of Alcoholism}
\description{
  Levels of expressed alpha synuclein mRNA in three groups of allele lengths of
  NACP-REP1.
}
\usage{alpha}
\format{
  A data frame with 97 observations on 2 variables.
  \describe{
    \item{\code{alength}}{
      allele length, a factor with levels \code{"short"}, \code{"intermediate"}
      and \code{"long"}.
    }
    \item{\code{elevel}}{
      expression levels of alpha synuclein mRNA.
    }
  }
}
\details{
  Various studies have linked alcohol dependence phenotypes to chromosome 4.
  One candidate gene is NACP (non-amyloid component of plaques), coding for
  alpha synuclein.  \bibcitet{coin::boenscheta2005} found longer
  alleles of NACP-REP1 in alcohol-dependent patients compared with healthy
  controls and reported that the allele lengths show some association with
  levels of expressed alpha synuclein mRNA.
}
\source{\bibshow{*}}
\examples{
## Boxplots
boxplot(elevel ~ alength, data = alpha)

## Asymptotic Kruskal-Wallis test
kruskal_test(elevel ~ alength, data = alpha)

## Asymptotic Kruskal-Wallis test using midpoint scores
kruskal_test(elevel ~ alength, data = alpha,
             scores = list(alength = c(2, 7, 11)))

## Asymptotic score-independent test
## Winell and Lindbaeck (2018)
(it <- independence_test(elevel ~ alength, data = alpha,
                         ytrafo = function(data)
                             trafo(data, numeric_trafo = rank_trafo), 
                         xtrafo = function(data)
                             trafo(data, factor_trafo = function(x)
                                 zheng_trafo(as.ordered(x)))))

## Extract the "best" set of scores
ss <- statistic(it, type = "standardized")
idx <- which(abs(ss) == max(abs(ss)), arr.ind = TRUE)
ss[idx[1], idx[2], drop = FALSE]
}
\keyword{datasets}
