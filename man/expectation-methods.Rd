\name{expectation-methods}
\docType{methods}
\alias{expectation}
\alias{expectation-methods}
\alias{expectation,IndependenceTest-method}
\alias{expectation,ScalarIndependenceTestStatistic-method}
\alias{expectation,MaxTypeIndependenceTestStatistic-method}
\alias{expectation,QuadTypeIndependenceTestStatistic-method}
\alias{covariance}
\alias{covariance-methods}
\alias{covariance,IndependenceTest-method}
\alias{covariance,ScalarIndependenceTestStatistic-method}
\alias{covariance,MaxTypeIndependenceTestStatistic-method}
\alias{covariance,QuadTypeIndependenceTestStatistic-method}
\title{ Extract the Expectation and Covariance of Linear Statistics }
\description{

  Extracts the conditional expectation and covariance 
  for linear statistics from objects inheriting from "IndependenceTest".

}
\usage{
expectation(object, ...)
covariance(object, ...)
}
\arguments{
  \item{object}{an object inheriting from class
    \code{\link{IndependenceTest-class}}.}
  \item{\dots}{further arguments (currently ignored).}
}
\section{Methods}{
  \describe{
    \item{expectation}{extracts the expectation of the linear 
      statistic of \code{object}.}
    \item{covariance}{extracts the covariance of the linear 
      statistic of \code{object}.}
  }
}
\examples{

df <- data.frame(y = gl(3, 2), x = gl(3, 2)[sample(1:6)])

### Cochran-Mantel-Haenzel Test
ct <- cmh_test(y ~ x, data = df)
 
### the linear statistic, i.e, the contingency table
l <- statistic(ct, type = "linear")
l

### expectation
El <- expectation(ct)
El

### covariance
Vl <- covariance(ct)
Vl

### the standardized contingency table (hard way)
(l - El) / matrix(sqrt(diag(Vl)), ncol = nrow(El))

### easy way
statistic(ct, type = "standardized")

}
\keyword{methods}
