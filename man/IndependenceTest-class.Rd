\name{IndependenceTest-class}
\docType{class}
\alias{IndependenceTest-class}
\alias{initialize,IndependenceTest-method}
\alias{ScalarIndependenceTest-class}
\alias{initialize,ScalarIndependenceTest-method}
\alias{ScalarIndependenceTestConfint-class}
\alias{initialize,ScalarIndependenceTestConfint-method}
\alias{MaxTypeIndependenceTest-class}
\alias{initialize,MayTypeIndependenceTest-method}
\alias{QuadTypeIndependenceTest-class}
\alias{initialize,QuadTypeIndependenceTest-method}

\title{Class "IndependenceTest"}
\description{  Objects of class "IndependenceTest" represent the results of
independence tests, including data, transformations, linear and test
statistics as well as their null distribution. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IndependenceTest", ...)}.
Normally, this class is only used internally.
}
\section{Slots}{
  \describe{
    \item{\code{distribution}:}{Object of class \code{"NullDistribution"}.}
    \item{\code{statistic}:}{Object of class \code{"IndependenceTestStatistic"}.}
    \item{\code{estimates}:}{Object of class \code{"list"}.}
    \item{\code{method}:}{Object of class \code{"character"}.}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "IndependenceTest")}: ...}
  }
}
\keyword{classes}
