\name{jobsatisfaction}
\docType{data}
\alias{jobsatisfaction}
\title{Income and Job Satisfaction}
\description{
  Income and job satisfaction by gender.
}
\usage{jobsatisfaction}
\format{
  A contingency table with 104 observations on 3 variables.
  \describe{
    \item{\code{Income}}{
      a factor with levels \code{"<5000"}, \code{"5000-15000"},
      \code{"15000-25000"} and \code{">25000"}.
    }
    \item{\code{Job.Satisfaction}}{
      a factor with levels \code{"Very Dissatisfied"},
      \code{"A Little Satisfied"}, \code{"Moderately Satisfied"} and
      \code{"Very Satisfied"}.
    }
    \item{\code{Gender}}{
      a factor with levels \code{"Female"} and \code{"Male"}.
    }
  }
}
\details{
  The data are given in Agresti (2002, p. 288, Tab. 7.8).
}
\source{
  Agresti, A.  (2002).  \emph{Categorical Data Analysis}, Second Edition.
  Hoboken, New Jersey: John Wiley & Sons.
}
\examples{
## Approximative (Monte Carlo) linear-by-linear association test
lbl_test(jobsatisfaction, distribution = approximate(B = 10000))
}
\keyword{datasets}
