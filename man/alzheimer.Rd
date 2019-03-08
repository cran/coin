\name{alzheimer}
\docType{data}
\alias{alzheimer}
\title{Smoking and Alzheimer's Disease}
\description{
  A case-control study of smoking and Alzheimer's disease.
}
\usage{alzheimer}
\format{
  A data frame with 538 observations on 3 variables.
  \describe{
    \item{\code{smoking}}{
      a factor with levels \code{"None"}, \code{"<10"}, \code{"10-20"} and
      \code{">20"} (cigarettes per day).
    }
    \item{\code{disease}}{
      a factor with levels \code{"Alzheimer"}, \code{"Other dementias"} and
      \code{"Other diagnoses"}.
    }
    \item{\code{gender}}{
      a factor with levels \code{"Female"} and \code{"Male"}.
    }
  }
}
\details{
  Subjects with Alzheimer's disease are compared to two different control groups
  with respect to smoking history.  The data are given in Salib and Hillier
  (1997, Tab. 4).
}
\source{
  Salib, E. and Hillier, V.  (1997).  A case-control study of smoking and
  Alzheimer's disease.  \emph{International Journal of Geriatric Psychiatry}
  \bold{12}(3), 295--300.
  \doi{10.1002/(SICI)1099-1166(199703)12:3<295::AID-GPS476>3.0.CO;2-3}
}
\references{
  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2006).  A Lego
  system for conditional inference.  \emph{The American Statistician}
  \bold{60}(3), 257--263.  \doi{10.1198/000313006X118430}
}
\examples{
## Spineplots
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:2, ncol = 2))
spineplot(disease ~ smoking, data = alzheimer,
          subset = gender == "Male", main = "Male")
spineplot(disease ~ smoking, data = alzheimer,
          subset = gender == "Female", main = "Female")
par(op) # reset

## Asymptotic Cochran-Mantel-Haenszel test
cmh_test(disease ~ smoking | gender, data = alzheimer)
}
\keyword{datasets}
