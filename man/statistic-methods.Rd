\name{statistic-methods}
\docType{methods}
\alias{statistic}
\alias{statistic-methods}
\alias{statistic,IndependenceLinearStatistic-method}
\alias{statistic,IndependenceTestStatistic-method}
\alias{statistic,IndependenceTest-method}
\title{Extraction of the Test Statistic and Linear Statistic}
\description{
  Methods for extraction of the test statistic and linear statistic.
}
\usage{
\S4method{statistic}{IndependenceLinearStatistic}(object, type = c("test", "linear", "centered", "standardized"), \dots)
\S4method{statistic}{IndependenceTestStatistic}(object, type = c("test", "linear", "centered", "standardized"), \dots)
\S4method{statistic}{IndependenceTest}(object, type = c("test", "linear", "centered", "standardized"), \dots)
}
\arguments{
  \item{object}{
    an object from which the test statistic or linear statistic can be
    extracted.
  }
  \item{type}{
    a character, the type of statistic: either \code{"test"} (default) for the
    test statistic, \code{"linear"} for the unstandardized linear statistic
    \code{"centered"} for the centered linear statistic or \code{"standardized"}
    for the standardized linear statistic.
  }
  \item{\dots}{
    further arguments (currently ignored).
  }
}
\details{
  The method \code{statistic} extracts the test statistic or the, possibly
  multivariate, linear statistic in its unstandardized, centered or
  standardized form.

  The test statistic (\code{type = "test"}) is returned by default.  The
  unstandardized, centered or standardized linear statistic is obtained by
  setting \code{type} to \code{"linear"}, \code{"centered"} or
  \code{"standardized"} respectively.
}
\value{
  The test statistic or the unstandardized, centered or standardized linear
  statistic extracted from \code{object}.  A numeric vector or matrix.
}
\examples{
## Example data
dta <- data.frame(
    y = gl(4, 5),
    x = gl(5, 4)
)

## Asymptotic Cochran-Mantel-Haenszel Test
ct <- cmh_test(y ~ x, data = dta)

## Test statistic
statistic(ct)

## The unstandardized linear statistic...
statistic(ct, type = "linear")

## ...is identical to the contingency table
xtabs(~ x + y, data = dta)

## The centered linear statistic...
statistic(ct, type = "centered")

## ...is identical to
statistic(ct, type = "linear") - expectation(ct)

## The standardized linear statistic, illustrating departures from the null
## hypothesis of independence...
statistic(ct, type = "standardized")

## ...is identical to
(statistic(ct, type = "linear") - expectation(ct)) / sqrt(variance(ct))
}
\keyword{methods}
