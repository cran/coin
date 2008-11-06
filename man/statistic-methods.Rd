\name{statistic-methods}
\docType{methods}
\alias{statistic}
\alias{statistic-methods}
\alias{statistic,IndependenceTest-method}
\alias{statistic,IndependenceLinearStatistic-method}
\alias{statistic,IndependenceTestStatistic-method}
\title{ Extract Test Statistics, Linear Statistics and Standardized
  Statistics }
\description{
  Extract the test statistic and, possibly multivariate, linear statistics
  in their raw and standardized form from objects of class
  \code{\link{IndependenceTest-class}}.
}
\usage{
statistic(object, type = c("test", "linear", "standardized"), ...)
}
\arguments{
  \item{object}{an object inheriting from class
    \code{\link{IndependenceTest-class}}.}
  \item{type}{either \code{test} for test statistic, \code{linear}
    for the unstandardized linear statistic or \code{standardized}
    for the standardized linear statistic.}
  \item{\dots}{further arguments (currently ignored).}
}
\section{Methods}{
  \describe{
    \item{statistic}{extracts the specified statistic from \code{object}.}
  }
}
\examples{

  df <- data.frame(y = gl(4, 5), x = gl(5, 4))

  ### Cochran-Mantel-Haenzel Test
  ct <- cmh_test(y ~ x, data = df)

  ### chi-squared statistic
  statistic(ct)

  ### the linear statistic, i.e, the contingency table
  statistic(ct, type = "linear")

  ### the same
  table(df$x, df$y)

  ### and the standardized contingency table for illustrating
  ### departures from the null hypothesis of independence of x and y
  statistic(ct, type = "standardized")

}
\keyword{methods}
