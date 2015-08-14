\name{IndependenceTestProblem-class}
\docType{class}
\alias{IndependenceTestProblem-class}
\alias{initialize,IndependenceTestProblem-method}
\title{Class \code{"IndependenceTestProblem"}}
\description{
  Objects of class \code{"IndependenceTestProblem"} represent the transformed
  and original data structures corresponding to an independence problem.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{     new("IndependenceTestProblem", object, xtrafo = trafo, ytrafo = trafo, \dots)}
  where \code{object} is an object of class
  \code{"\linkS4class{IndependenceProblem}"}, \code{xtrafo} is the regression
  function \eqn{g(\mathbf{X})}{g(X)} and \code{ytrafo} is the influence function
  \eqn{h(\mathbf{Y})}{h(Y)}.
}
\section{Slots}{
  \describe{
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
  Class \code{"\linkS4class{IndependenceProblem}"}, directly.
}
\section{Known Subclasses}{
  Class \code{"\linkS4class{IndependenceLinearStatistic}"}, directly. \cr
  Class \code{"\linkS4class{IndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceLinearStatistic}"}, distance 2. \cr
  Class \code{"\linkS4class{MaxTypeIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 3. \cr
  Class \code{"\linkS4class{QuadTypeIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 3. \cr
  Class \code{"\linkS4class{ScalarIndependenceTestStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 3.
}
\section{Methods}{
  \describe{
    \item{initialize}{
      \code{signature(.Object = "IndependenceTestProblem")}: See the
      documentation for \code{\link[methods:new]{initialize}} (in
      package \pkg{methods}) for details.
    }
  }
}
\keyword{classes}
