\name{neuropathy}
\docType{data}
\alias{neuropathy}
\title{Acute Painful Diabetic Neuropathy}
\description{
  The logarithm of the ratio of pain scores measured at baseline and after four
  weeks in a control group and a treatment group.
}
\usage{neuropathy}
\format{
  A data frame with 58 observations on 2 variables.
  \describe{
    \item{\code{pain}}{
      pain scores: ln(baseline / final).
    }
    \item{\code{group}}{
      a factor with levels \code{"control"} and \code{"treat"}.
    }
  }
}
\details{
  Data from Conover and Salsburg (1988, Tab. 1).
}
\source{
  Conover, W. J. and Salsburg, D. S.  (1988).  Locally most powerful tests for
  detecting treatment effects when only a subset of patients can be expected to
  \dQuote{respond} to treatment.  \emph{Biometrics} \bold{44}(1), 189--196.
}
\examples{
## Conover and Salsburg (1988, Tab. 2)

## One-sided approximative Fisher-Pitman test
oneway_test(pain ~ group, data = neuropathy,
            alternative = "less",
            distribution = approximate(B = 10000))

## One-sided approximative Wilcoxon-Mann-Whitney test
wilcox_test(pain ~ group, data = neuropathy,
            alternative = "less",
            distribution = approximate(B = 10000))

## One-sided approximative Conover-Salsburg test
oneway_test(pain ~ group, data = neuropathy,
            alternative = "less",
            distribution = approximate(B = 10000),
            ytrafo = function(data)
                trafo(data, numeric_trafo = consal_trafo))

## One-sided approximative maximum test for a range of 'a' values
it <- independence_test(pain ~ group, data = neuropathy,
                        alternative = "less",
                        distribution = approximate(B = 10000),
                        ytrafo = function(data)
                            trafo(data, numeric_trafo = function(y)
                                consal_trafo(y, a = 2:7)))
pvalue(it, method = "single-step")
}
\keyword{datasets}
