\name{rotarod}
\alias{rotarod}
\docType{data}
\title{ Rotating Rats Data}
\description{
  The endurance time of 24 rats in two groups in a rotating cylinder.
}
\usage{data("rotarod")}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{time}{the endurance time}
    \item{group}{a factor with levels \code{control} and \code{treatment}.}}
}
\details{
The 24 rats received a fixed oral dose of a centrally acting muscle relaxant
(\code{treatment}) or a saline solvent (\code{control}). They were placed on a rotating
cylinder and the length of time each rat remains on the cylinder is
measured, up to a maximum of 300 seconds.
The rats were randomly assigned to the control and treatment group.

Note that the empirical variance in the control group is 0 and that the
group medians are identical.

This dataset serves as the basis of an comparision of the results of the 
Wilcoxon-Mann-Whitney test computed by 11 statistical packages in Bergmann
et al. (2000). The exact conditional p-value is 0.0373 (two-sided) and
0.0186 (one-sided). The asymptotic two-sided 
p-value (corrected for ties) is reported as 0.0147.
}
\source{

    Reinhard Bergmann, John Ludbrook & Will P. J. M. Spooren (2000),
    Different outcomes of the Wilcoxon-Mann-Whitney test from different 
    statistics packages. \emph{The American Statistician}
    \bold{54}(1), 72--77.

}
\examples{

data("rotarod", package = "coin")

### Wilcoxon-Mann-Whitney Rank Sum Test

### one-sided exact (0.0186)
wilcox_test(time ~ group, data = rotarod, 
    alternative = "greater", distribution = "exact")
### two-sided exact (0.0373)
wilcox_test(time ~ group, data = rotarod, distribution = "exact")
### two-sided asymptotical (0.0147)
wilcox_test(time ~ group, data = rotarod)

}
\keyword{datasets}
