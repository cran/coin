\name{alzheimer}
\alias{alzheimer}
\docType{data}
\title{ Smoking and Alzheimer's Disease }
\description{
  A case-control study of smoking and Alzheimer's disease.
}
\usage{data("alzheimer")}
\format{
 A tree-dimensional contingency table with dimensions
 \description{
     \item{smoking}{a factor at levels \code{None}, \code{<10}, \code{10-20}
                    and \code{>20} (cigarettes per day).}
     \item{disease}{a factor at levels \code{Alzheimer's}, \code{Other
                    dementias} and \code{Other diagnoses}.}
     \item{gender}{a factor at levels \code{Female} and \code{Male}.}
 }
}
\details{
    198 cases of Alzheimer's disease are compared to a control group with
    respect to smoking history. The data are published in Table 4 (Salib \&
    Hillier, 1997).
}
\source{

    Emad Salib \& Valerie Hillier (1997), A case-control study of smoking
    and Alzheimer's disease. \emph{International Journal of Geriatric
    Psychiatry}, \bold{12}, 295--300.

}
\examples{
data("alzheimer")

### Cochran-Mantel-Haenszel test
cmh_test(alzheimer)

### Linear-by-Linear Association test
cmh_test(alzheimer, scores = list(smoking = c(0, 5, 15, 25)))
statistic(cmh_test(alzheimer, scores = list(smoking = c(0, 5, 15, 25))),
          "standardized")

}
\keyword{datasets}
