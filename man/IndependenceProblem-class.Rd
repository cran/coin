\name{IndependenceProblem-class}
\docType{class}
\alias{IndependenceProblem-class}
\alias{initialize,IndependenceProblem-method}
\title{Class \code{"IndependenceProblem"}}
\description{
  Objects of class \code{"IndependenceProblem"} represent the data structure
  corresponding to an independence problem.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{     new("IndependenceProblem", x, y, block = NULL, weights = NULL, \dots)}
  where \code{x} and \code{y} are data frames containing the variables
  \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y} respectively, \code{block} is an
  optional factor representing the block structure \eqn{b} and \code{weights} is
  an optional integer vector corresponding to the case weights \eqn{w}.
}
\section{Slots}{
  \describe{
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
\section{Known Subclasses}{
  Class \code{"\linkS4class{IndependenceTestProblem}"}, directly. \cr
  Class \code{"\linkS4class{SymmetryProblem}"}, directly. \cr
  Class \code{"\linkS4class{IndependenceLinearStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestProblem}"}, distance 2. \cr
  Class \code{"\linkS4class{IndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestProblem}"}, distance 3. \cr
  Class \code{"\linkS4class{MaxTypeIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestProblem}"}, distance 4. \cr
  Class \code{"\linkS4class{QuadTypeIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestProblem}"}, distance 4. \cr
  Class \code{"\linkS4class{ScalarIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestProblem}"}, distance 4.
}
\section{Methods}{
  \describe{
    \item{initialize}{
      \code{signature(.Object = "IndependenceProblem")}: See the documentation
      for \code{\link[methods:new]{initialize}} (in package \pkg{methods}) for
      details.
    }
  }
}
\keyword{classes}
