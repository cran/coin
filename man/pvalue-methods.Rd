\name{pvalue-methods}
\docType{methods}
\alias{pvalue}
\alias{pvalue-methods}
\alias{pvalue,NullDistribution-method}
\alias{pvalue,IndependenceTest-method}
\alias{pvalue,ScalarIndependenceTest-method}
\alias{pvalue,MaxTypeIndependenceTest-method}
\alias{pvalue,QuadTypeIndependenceTest-method}
\title{ Extract P-Values }
\description{

  Extracts the p-value from objects representing null distributions of 
  independence tests.

}
\usage{
pvalue(object, ...)
}
\arguments{
  \item{object}{an object inheriting from class
    \code{\link{IndependenceTest-class}}.}
  \item{\dots}{further arguments (currently ignored).}
}
\section{Methods}{
  \describe{
    \item{pvalue}{extracts the p-value from the specified object.}
  }
}
\examples{

### artificial 2-sample problem
df <- data.frame(y = rnorm(20), x = gl(2, 10))
 
### Ansari-Bradley test
at <- ansari_test(y ~ x, data = df, distribution = "exact")

at

pvalue(at)

}
\keyword{methods}
\keyword{htest}
