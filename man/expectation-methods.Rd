\name{expectation-methods}
\docType{methods}
\alias{expectation}
\alias{expectation-methods}
\alias{expectation,IndependenceLinearStatistic-method}
\alias{expectation,IndependenceTest-method}
\alias{variance}
\alias{variance-methods}
\alias{variance,Variance-method}
\alias{variance,CovarianceMatrix-method}
\alias{variance,IndependenceLinearStatistic-method}
\alias{variance,IndependenceTest-method}
\alias{covariance}
\alias{covariance-methods}
\alias{covariance,CovarianceMatrix-method}
\alias{covariance,IndependenceLinearStatistic-method}
\alias{covariance,IndependenceTest-method}
\title{Extraction of the Expectation, Variance and Covariance of the Linear
  Statistic}
\description{
  Methods for extraction of the expectation, variance and covariance of the
  linear statistic.
}
\usage{
\S4method{expectation}{IndependenceLinearStatistic}(object, \dots)
\S4method{expectation}{IndependenceTest}(object, \dots)

\S4method{variance}{Variance}(object, \dots)
\S4method{variance}{CovarianceMatrix}(object, \dots)
\S4method{variance}{IndependenceLinearStatistic}(object, \dots)
\S4method{variance}{IndependenceTest}(object, \dots)

\S4method{covariance}{CovarianceMatrix}(object, \dots)
\S4method{covariance}{IndependenceLinearStatistic}(object, \dots)
\S4method{covariance}{IndependenceTest}(object, \dots)
}
\arguments{
  \item{object}{
    an object from which the expectation, variance or covariance of the linear
    statistic can be extracted.
  }
  \item{\dots}{
    further arguments (currently ignored).
  }
}
\details{
  The methods \code{expectation}, \code{variance} and \code{covariance} extract
  the expectation, variance and covariance, respectively, of the linear
  statistic.
}
\value{
  The expectation, variance or covariance of the linear statistic extracted from
  \code{object}.  A numeric vector or matrix.
}
\examples{\dontshow{suppressWarnings(RNGversion("3.5.2")); set.seed(711109)}
## Example data
dta <- data.frame(
    y = gl(3, 2),
    x = sample(gl(3, 2))
)

## Asymptotic Cochran-Mantel-Haenszel Test
ct <- cmh_test(y ~ x, data = dta)

## The linear statistic, i.e., the contingency table...
(l <- statistic(ct, type = "linear"))

## ...and its expectation...
(El <- expectation(ct))

## ...and covariance
(Vl <- covariance(ct))

## The standardized contingency table...
(l - El) / sqrt(variance(ct))

## ...is identical to the standardized linear statistic
statistic(ct, type = "standardized")
}
\keyword{methods}
