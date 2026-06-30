\name{asat}
\docType{data}
\alias{asat}
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
  compound.  The data were originally given in \bibcitet{coin::Hothorn1992} and later
  reproduced by \bibcitet{coin::Hauschke_1999}, see also
  \bibcitet{coin::pflueger_2002}.
}
\references{\bibshow{*}}
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
