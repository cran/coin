\name{asat}
\alias{asat}
\docType{data}
\title{ Toxicological Study on Female Wistar Rats }
\description{
  ASAT-values for a new compound and a control group of 34 female Wistar
rats.
}
\usage{data("asat")}
\format{
  A data frame with 34 observations on the following 2 variables.
  \describe{
    \item{asat}{the ASAT-values (a liver enzyme)}
    \item{group}{a factor with levels \code{Compound} and \code{Control}.}
  }
}
\details{
  The aim of this toxicological study is the proof of safety for the new
compound. The data are originally given in Hothorn (1992) and reproduced in
Hauschke et al. (1999).
}
\source{

  Ludwig A. Hothorn (1992), Biometrische Analyse toxikologischer Untersuchungen.
  In: J. Adams (ed.): \emph{Statistisches Know how in der medizinischen 
  Forschung.} Ullstein-Mosby, Berlin, 475--590.

}
\references{

  Dieter Hauschke, Meinhard Kieser \& Ludwig A. Hothorn (1999).
  Proof of safety in
  toxicology based on the ratio of two means for normally distributed data,
  \emph{Biometrical Journal} \bold{41}(3), 295--304.

  Rafael Pfl{\"u}ger \& Torsten Hothorn (2002).
  Assessing equivalence tests with respect to their expected
  $p$-Value, \emph{Biometrical Journal} \bold{44}(8), 1002--1027.

}
\examples{

  ### proof-of-safety based on ratio of medians
  pos <- wilcox_test(I(log(asat)) ~ group, data = asat, alternative = "less", 
                     conf.int = TRUE, distribution = "exact")

  ### one-sided confidence set. Safety cannot be concluded since the effect of
  ### the compound exceeds 20\% of the control median
  exp(confint(pos)$conf.int)

}
\keyword{datasets}
