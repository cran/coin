\name{rotarod}
\docType{data}
\alias{rotarod}
\title{Rotating Rats}
\description{
  The endurance time of 24 rats in two groups on a rotating cylinder.
}
\usage{rotarod}
\format{
  A data frame with 24 observations on 2 variables.
  \describe{
    \item{\code{time}}{
      endurance time (seconds).
    }
    \item{\code{group}}{
      a factor with levels \code{"control"} and \code{"treatment"}.
    }
  }
}
\details{
  The rats were randomly assigned to receive a fixed oral dose of a centrally
  acting muscle relaxant (\code{"treatment"}) or a saline solvent
  (\code{"control"}).  The animals were placed on a rotating cylinder and the
  endurance time of each rat, i.e., the length of time each rat remained on the
  cylinder, was measured up to a maximum of 300 seconds.

  This dataset is the basis of a comparison of 11 different software
  implementations of the Wilcoxon-Mann-Whitney test presented in Bergmann,
  Ludbrook and Spooren (2000).
}
\note{
  The empirical variance in the control group is 0 and the group medians are
  identical.  The exact conditional \eqn{p}-values are 0.0373 (two-sided) and
  0.0186 (one-sided).  The asymptotic two-sided \eqn{p}-value (corrected for
  ties) is 0.0147.
}
\source{
  Bergmann, R., Ludbrook, J. and Spooren, W. P. J. M.  (2000).  Different
  outcomes of the Wilcoxon-Mann-Whitney test from different statistics packages.
  \emph{The American Statistician} \bold{54}(1), 72--77.
}
\examples{
## One-sided exact Wilcoxon-Mann-Whitney test (p = 0.0186)
wilcox_test(time ~ group, data = rotarod, distribution = "exact",
            alternative = "greater")

## Two-sided exact Wilcoxon-Mann-Whitney test (p = 0.0373)
wilcox_test(time ~ group, data = rotarod, distribution = "exact")

## Two-sided asymptotic Wilcoxon-Mann-Whitney test (p = 0.0147)
wilcox_test(time ~ group, data = rotarod)
}
\keyword{datasets}
