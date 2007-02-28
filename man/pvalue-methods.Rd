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
  \item{\dots}{additional arguments: \code{method}, 
      a character specifying the type of 
      adjustment (\code{global}, \code{single-step}, \code{step-down} or \code{discrete})
      should be used. The default is \code{global}.}
}
\section{Methods}{
  \describe{
    \item{pvalue}{extracts the p-value from the specified object.}
  }
}
\details{
  
  Univariate p-values for maximum-type statistics come with associated 99\%
  confidence interval when resampling was used to determine the null
  distribution (which may be the case even when \code{distribution =
  "asypmtotic"} was used). 

  By default, a global p-value is returned. When \code{method =
  "single-step"}, adjusted p-values are obtained from a 
  single-step max-T procedure
  (Westfall & Young, 1993, algorithm 2.5 and formula 2.8). Note that the
  minimum of the adjusted p-values always controls the familywise error 
  rate (FWER) but the maximum type I error, i.e. the error for
  each of the individual tests, is only controlled when the subset
  pivotality condition holds.

  When \code{method = "step-down"} the free step-down resampling method
  (algorithm 2.8 and formula 2.8 in Westfall & Young, 1993) is used, the above
  comments apply as well.

  With \code{method = "discrete"}, the Bonferroni adjustment as suggested by
  Westfall & Wolfinger (1997) with improvements for highly discrete
  permutation distributions is available, however, without taking
  correlations between the test statistics into account. Here, the p-values are
  valid even without assuming subset pivotality.

}
\references{

    Peter H. Westfall \& S. Stanley Young (1993).
    \emph{Resampling-based Multiple Testing}.
    New York: John Wiley & Sons.

    Peter H. Westfall \& Russell D. Wolfinger (1997).
    Multiple tests with discrete distributions. 
    \emph{The American Statistician} \bold{51}, 3--8.

}
\examples{

  ### artificial 2-sample problem
  df <- data.frame(y = rnorm(20), x = gl(2, 10))
 
  ### Ansari-Bradley test
  at <- ansari_test(y ~ x, data = df, distribution = "exact")
  at

  pvalue(at)

  ### bivariate 2-sample problem
  df <- data.frame(y1 = rnorm(20) + c(rep(0, 10), rep(1, 10)), 
                   y2 = rnorm(20), 
                   x = gl(2, 10))

  it <- independence_test(y1 + y2 ~ x, data = df, 
                          distribution = approximate(B = 9999))
  pvalue(it, method = "single-step")
  pvalue(it, method = "step-down")
  pvalue(it, method = "discrete")

}
\keyword{methods}
\keyword{htest}
