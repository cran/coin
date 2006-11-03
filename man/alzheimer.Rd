\name{alzheimer}
\alias{alzheimer}
\docType{data}
\title{ Smoking and Alzheimer's Disease }
\description{
  A case-control study of smoking and Alzheimer's disease.
}
\usage{data("alzheimer")}
\format{
  A data frame with 538 observations on the following 3 variables.
 \describe{
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

    Emad Salib \& Valerie Hillier (1997). A case-control study of smoking
    and Alzheimer's disease, \emph{International Journal of Geriatric
    Psychiatry} \bold{12}, 295--300.

    Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel \& Achim Zeileis (2006).
    A Lego system for conditional inference, \emph{The American Statistician},
    \bold{60}(3), 257--263.


}
\examples{

  data("alzheimer", package = "coin")

  ### spineplots
  layout(matrix(1:2, ncol = 2))
  spineplot(disease ~ smoking, data = alzheimer, subset = gender == "Male",
            main = "Male")
  spineplot(disease ~ smoking, data = alzheimer, subset = gender == "Female",
            main = "Female")

  ### Cochran-Mantel-Haenszel test
  cmh_test(disease ~ smoking | gender, data = alzheimer)

}
\keyword{datasets}
