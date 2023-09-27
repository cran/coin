\name{NullDistribution-methods}
\docType{methods}
\alias{AsymptNullDistribution}
\alias{AsymptNullDistribution-methods}
\alias{AsymptNullDistribution,MaxTypeIndependenceTestStatistic-method}
\alias{AsymptNullDistribution,QuadTypeIndependenceTestStatistic-method}
\alias{AsymptNullDistribution,ScalarIndependenceTestStatistic-method}
\alias{ApproxNullDistribution}
\alias{ApproxNullDistribution-methods}
\alias{ApproxNullDistribution,MaxTypeIndependenceTestStatistic-method}
\alias{ApproxNullDistribution,QuadTypeIndependenceTestStatistic-method}
\alias{ApproxNullDistribution,ScalarIndependenceTestStatistic-method}
\alias{ExactNullDistribution}
\alias{ExactNullDistribution-methods}
\alias{ExactNullDistribution,QuadTypeIndependenceTestStatistic-method}
\alias{ExactNullDistribution,ScalarIndependenceTestStatistic-method}
\title{Computation of the Reference Distribution}
\description{
  Methods for computation of the asymptotic, approximative (Monte Carlo) and
  exact reference distribution.
}
\usage{
\S4method{AsymptNullDistribution}{MaxTypeIndependenceTestStatistic}(object, ...)
\S4method{AsymptNullDistribution}{QuadTypeIndependenceTestStatistic}(object, ...)
\S4method{AsymptNullDistribution}{ScalarIndependenceTestStatistic}(object, ...)

\S4method{ApproxNullDistribution}{MaxTypeIndependenceTestStatistic}(object, nresample = 10000L, B, ...)
\S4method{ApproxNullDistribution}{QuadTypeIndependenceTestStatistic}(object, nresample = 10000L, B, ...)
\S4method{ApproxNullDistribution}{ScalarIndependenceTestStatistic}(object, nresample = 10000L, B, ...)

\S4method{ExactNullDistribution}{QuadTypeIndependenceTestStatistic}(object, algorithm = c("auto", "shift", "split-up"), ...)
\S4method{ExactNullDistribution}{ScalarIndependenceTestStatistic}(object, algorithm = c("auto", "shift", "split-up"), ...)
}
\arguments{
  \item{object}{
    an object from which the asymptotic, approximative (Monte Carlo) or exact
    reference distribution can be computed.
  }
  \item{nresample}{
    a positive integer, the number of Monte Carlo replicates used for the
    computation of the approximative reference distribution.  Defaults to
    \code{10000L}.
  }
  \item{B}{
    deprecated, use \code{nresample} instead.
  }
  \item{algorithm}{
    a character, the algorithm used for the computation of the exact reference
    distribution: either \code{"auto"} (default), \code{"shift"} or
    \code{"split-up"}.
  }
  \item{...}{
    further arguments to be passed to or from methods.
  }
}
\details{
  The methods \code{AsymptNullDistribution}, \code{ApproxNullDistribution} and
  \code{ExactNullDistribution} compute the asymptotic, approximative (Monte
  Carlo) and exact reference distribution, respectively.
}
\note{
  In versions of \pkg{coin} prior to 1.3-0, the number of Monte Carlo replicates
  in \code{ApproxNullDistribution()} was specified using the now deprecated
  \code{B} argument.  \strong{This will be made defunct and removed in a future
  release.}  It has been replaced by the \code{nresample} argument (for
  conformity with the \pkg{libcoin}, \pkg{party} and \pkg{partykit} packages).
}
\value{
  An object of class \code{"\linkS4class{AsymptNullDistribution}"},
  \code{"\linkS4class{ApproxNullDistribution}"} or
  \code{"\linkS4class{ExactNullDistribution}"}.
}
\keyword{methods}
