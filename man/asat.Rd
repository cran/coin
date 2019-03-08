\name{asat}
\docType{data}
\alias{asat}
\encoding{UTF-8}
\title{Toxicological Study on Female Wistar Rats}
\description{
  Measurements of the liver enzyme aspartate aminotransferase (ASAT) for a new
  compound and a control group of 34 female Wistar rats.
}
\usage{asat}
\format{
  A data frame with 34 observations on 2 variables.
  \describe{
    \item{\code{asat}}{
      ASAT values.
    }
    \item{\code{group}}{
      a factor with levels \code{"Compound"} and \code{"Control"}.
    }
  }
}
\details{
  The aim of this toxicological study is the proof of safety for the new
  compound.  The data were originally given in Hothorn (1992) and later
  reproduced by Hauschke, Kieser and Hothorn (1999).
}
\source{
  Hauschke, D., Kieser, M. and Hothorn, L. A.  (1999).  Proof of safety in
  toxicology based on the ratio of two means for normally distributed data.
  \emph{Biometrical Journal} \bold{41}(3), 295--304.
  \doi{10.1002/(SICI)1521-4036(199906)41:3<295::AID-BIMJ295>3.0.CO;2-2}

  Hothorn, L. A.  (1992).  Biometrische analyse toxikologischer untersuchungen.
  In J. Adam (Ed.), \emph{Statistisches Know-How in der Medizinischen
    Forschung}, pp. 475--590.  Berlin: Ullstein Mosby.
}
\references{
  \enc{Pflüger}{Pflueger}, R. and Hothorn, T.  (2002).  Assessing equivalence
  tests with respect to their expected \eqn{p}-value.  \emph{Biometrical
    Journal} \bold{44}(8), 1015--1027.  \doi{10.1002/bimj.200290001}
}
\examples{
## Proof-of-safety based on ratio of medians (Pflueger and Hothorn, 2002)
## One-sided exact Wilcoxon-Mann-Whitney test
wt <- wilcox_test(I(log(asat)) ~ group, data = asat,
                  distribution = "exact", alternative = "less",
                  conf.int = TRUE)

## One-sided confidence set
## Note: Safety cannot be concluded since the effect of the compound
##       exceeds 20 \% of the control median
exp(confint(wt)$conf.int)
}
\keyword{datasets}
