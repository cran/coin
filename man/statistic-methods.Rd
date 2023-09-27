\name{statistic-methods}
\docType{methods}
\alias{statistic}
\alias{statistic-methods}
\alias{statistic,IndependenceLinearStatistic-method}
\alias{statistic,IndependenceTestStatistic-method}
\alias{statistic,IndependenceTest-method}
\title{Extraction of the Test Statistic and the Linear Statistic}
\description{
  Methods for extraction of the test statistic and the linear statistic.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\S4method{statistic}{IndependenceLinearStatistic}(object, type = c("test", "linear", "centered", "standardized"),
          partial = FALSE, ...)
\S4method{statistic}{IndependenceTestStatistic}(object, type = c("test", "linear", "centered", "standardized"),
          partial = FALSE, ...)
\S4method{statistic}{IndependenceTest}(object, type = c("test", "linear", "centered", "standardized"),
          partial = FALSE, ...)
}
\arguments{
  \item{object}{
    an object from which the test statistic or the linear statistic can be
    extracted.
  }
  \item{type}{
    a character string indicating the type of statistic: either \code{"test"}
    (default) for the test statistic, \code{"linear"} for the unstandardized
    linear statistic, \code{"centered"} for the centered linear statistic or
    \code{"standardized"} for the standardized linear statistic.
  }
  \item{partial}{
    a logical indicating that the partial linear statistic for each block should
    be extracted.  Defaults to \code{FALSE}.
  }
  \item{...}{
    further arguments (currently ignored).
  }
}
\details{
  The method \code{statistic} extracts the univariate test statistic or the,
  possibly multivariate, linear statistic in its unstandardized, centered or
  standardized form.

  The test statistic (\code{type = "test"}) is returned by default.  The
  unstandardized, centered or standardized linear statistic is obtained by
  setting \code{type} to \code{"linear"}, \code{"centered"} or
  \code{"standardized"}, respectively.  For tests of conditional independence
  within blocks, the partial linear statistic for each block is obtained by
  setting \code{partial = TRUE}.
}
\value{
  The test statistic or the unstandardized, centered or standardized linear
  statistic extracted from \code{object}.  A numeric vector, matrix or array.
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
