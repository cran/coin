\name{alpha}
\docType{data}
\alias{alpha}
\encoding{UTF-8}
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
  alpha synuclein.  \enc{Bönsch}{Boensch} \emph{et al.} (2005) found longer
  alleles of NACP-REP1 in alcohol-dependent patients compared with healthy
  controls and reported that the allele lengths show some association with
  levels of expressed alpha synuclein mRNA.
}
\source{
  \enc{Bönsch}{Boensch}, D., Lederer, T., Reulbach, U., Hothorn, T., Kornhuber,
  J. and Bleich, S.  (2005).  Joint analysis of the \emph{NACP}-REP1 marker
  within the alpha synuclein gene concludes association with alcohol dependence.
  \emph{Human Molecular Genetics} \bold{14}(7), 967--971.
  \doi{10.1093/hmg/ddi090}
}
\references{
  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2006).  A Lego
  system for conditional inference.  \emph{The American Statistician}
  \bold{60}(3), 257--263.  \doi{10.1198/000313006X118430}
}
\examples{
## Boxplots
boxplot(elevel ~ alength, data = alpha)

## Asymptotic Kruskal-Wallis test
kruskal_test(elevel ~ alength, data = alpha)
}
\keyword{datasets}
