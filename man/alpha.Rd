\name{alpha}
\alias{alpha}
\docType{data}
\title{ Genetic Components of Alcoholism }
\description{
  Levels of expressed alpha synuclein mRNA in three groups of allele lengths
  of NACP-REP1.
}
\usage{data("alpha")}
\format{
  A data frame with 97 observations on the following 2 variables.
  \describe{
    \item{\code{alength}}{allele length, a factor with levels \code{short},
                         \code{intermediate} and \code{long}.}
    \item{\code{elevel}}{expression levels of alpha synuclein mRNA.}
  }
}
\details{

  Various studies have linked alcohol dependence phenotypes to chromosome 4.
  One candidate gene is NACP (non-amyloid component of plaques),
  coding for alpha synuclein. B\"onsch et al. (2005) found longer alleles of
  NACP-REP1 in alcohol-dependent patients compared with healthy
  controls and report that the allele lengths show some
  association with levels of expressed alpha synuclein mRNA.

}
\source{

  Dominikus B\"onsch, Thomas Lederer, Udo Reulbach,  Torsten Hothorn, 
  Johannes Kornhuber \& Stefan Bleich (2005).
  Joint Analysis of the NACP-REP1 Marker Within the Alpha
  Synuclein Gene Concludes Association with Alcohol Dependence.
  \emph{Human Molecular Genetics}, \bold{14}(7), 967--971.

  Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel \& Achim Zeileis (2006).
  A Lego system for conditional inference, \emph{The American Statistician},
  \bold{60}(3), 257--263.

}
\examples{

  boxplot(elevel ~ alength, data = alpha)
  kruskal_test(elevel ~ alength, data = alpha)

}
\keyword{datasets}
