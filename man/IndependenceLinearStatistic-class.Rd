\name{IndependenceLinearStatistic-class}
\docType{class}
\alias{IndependenceLinearStatistic-class}
\alias{initialize,IndependenceLinearStatistic-method}
\title{Class \code{"IndependenceLinearStatistic"}}
\description{
  Objects of class \code{"IndependenceLinearStatistic"} represent the linear
  statistic and the transformed and original data structures corresponding to an
  independence problem.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{     new("IndependenceLinearStatistic", object, varonly = FALSE, \dots)}
  where \code{object} is an object of class
  \code{"\linkS4class{IndependenceTestProblem}"}, \code{varonly} is a logical
  indicating that slot \code{covariance} (see \sQuote{Slots}) should only
  contain the diagonal elements of the covariance matrix.
}
\section{Slots}{
  \describe{
    \item{\code{linearstatistic}:}{
      Object of class \code{"numeric"}.  The linear statistic.
    }
    \item{\code{expectation}:}{
      Object of class \code{"numeric"}.  The expectation of the linear
      statistic.
    }
    \item{\code{covariance}:}{
      Object of class \code{"\linkS4class{VarCovar}"}.  The covariance or
      variance of the linear statistic.
    }
    \item{\code{xtrans}:}{
      Object of class \code{"matrix"}.  The transformed \code{x}.
    }
    \item{\code{ytrans}:}{
      Object of class \code{"matrix"}.  The transformed \code{y}.
    }
    \item{\code{xtrafo}:}{
      Object of class \code{"function"}.  The regression function for \code{x}.
    }
    \item{\code{ytrafo}:}{
      Object of class \code{"function"}.  The influence function for \code{y}.
    }
    \item{\code{x}:}{
      Object of class \code{"data.frame"}.  The variables \code{x}.
    }
    \item{\code{y}:}{
      Object of class \code{"data.frame"}.  The variables \code{y}.
    }
    \item{\code{block}:}{
      Object of class \code{"factor"}.  The block structure.
    }
    \item{\code{weights}:}{
      Object of class \code{"numeric"}.  The case weights.
    }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{IndependenceTestProblem}"}, directly. \cr
  Class \code{"\linkS4class{IndependenceProblem}"}, by class
  \code{"\linkS4class{IndependenceTestProblem}"}, distance 2.
}
\section{Known Subclasses}{
  Class \code{"\linkS4class{IndependenceTestStatistic}"}, directly.\cr
  Class \code{"\linkS4class{MaxTypeIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 2. \cr
  Class \code{"\linkS4class{QuadTypeIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 2. \cr
  Class \code{"\linkS4class{ScalarIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 2.
}
\section{Methods}{
  \describe{
    \item{covariance}{
      \code{signature(object = "IndependenceLinearStatistic")}: See the
      documentation for \code{\link{covariance}} for details.
    }
    \item{expectation}{
      \code{signature(object = "IndependenceLinearStatistic")}: See the
      documentation for \code{\link{expectation}} for details.
    }
    \item{initialize}{
      \code{signature(.Object = "IndependenceLinearStatistic")}: See the
      documentation for \code{\link[methods:new]{initialize}} (in package
      \pkg{methods}) for details.
    }
    \item{statistic}{
      \code{signature(object = "IndependenceLinearStatistic")}: See the
      documentation for \code{\link{statistic}} for details.
    }
    \item{variance}{
      \code{signature(object = "IndependenceLinearStatistic")}: See the
      documentation for \code{\link{variance}} for details.
    }
  }
}
\keyword{classes}
