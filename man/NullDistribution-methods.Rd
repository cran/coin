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
\S4method{AsymptNullDistribution}{MaxTypeIndependenceTestStatistic}(object, \dots)
\S4method{AsymptNullDistribution}{QuadTypeIndependenceTestStatistic}(object, \dots)
\S4method{AsymptNullDistribution}{ScalarIndependenceTestStatistic}(object, \dots)

\S4method{ApproxNullDistribution}{MaxTypeIndependenceTestStatistic}(object, B = 10000, \dots)
\S4method{ApproxNullDistribution}{QuadTypeIndependenceTestStatistic}(object, B = 10000, \dots)
\S4method{ApproxNullDistribution}{ScalarIndependenceTestStatistic}(object, B = 10000, \dots)

\S4method{ExactNullDistribution}{QuadTypeIndependenceTestStatistic}(object, algorithm = c("auto", "shift", "split-up"), \dots)
\S4method{ExactNullDistribution}{ScalarIndependenceTestStatistic}(object, algorithm = c("auto", "shift", "split-up"), \dots)
}
\arguments{
  \item{object}{
    an object from which the asymptotic, approximative (Monte Carlo) or exact
    reference distribution can be computed.
  }
  \item{B}{
    a positive integer, the number of Monte Carlo replicates used for the
    computation of the approximative reference distribution.  Defaults to
    \code{10000}.
  }
  \item{algorithm}{
    a character, the algorithm used for the computation of the exact reference
    distribution: either \code{"auto"} (default), \code{"shift"} or
    \code{"split-up"}.
  }
  \item{\dots}{
    further arguments to be passed to or from methods.
  }
}
\details{
  The methods \code{AsymptNullDistribution}, \code{ApproxNullDistribution} and
  \code{ExactNullDistribution} compute the asymptotic, approximative (Monte
  Carlo) and exact reference distribution respectively.
}
\value{
  An object of class \code{"\linkS4class{AsymptNullDistribution}"},
  \code{"\linkS4class{ApproxNullDistribution}"} or
  \code{"\linkS4class{ExactNullDistribution}"}.
}
\keyword{methods}