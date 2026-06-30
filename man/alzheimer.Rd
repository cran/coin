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
  with respect to smoking history.  The data are given in
  \bibcitet{|salibhillier1997|Tab. 4}.
}
\source{\bibshow{*}}
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
