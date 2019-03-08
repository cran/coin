\name{SymmetryProblem-class}
\docType{class}
\alias{SymmetryProblem-class}
\alias{initialize,SymmetryProblem-method}
\title{Class \code{"SymmetryProblem"}}
\description{
  Objects of class \code{"SymmetryProblem"} represent the data structure
  corresponding to a symmetry problem.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{     new("SymmetryProblem", x, y, block = NULL, weights = NULL, \dots)}
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
\section{Extends}{
  Class \code{"\linkS4class{IndependenceProblem}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{
      \code{signature(.Object = "SymmetryProblem")}: See the documentation for
      \code{\link[methods:new]{initialize}} (in package \pkg{methods}) for
      details.
    }
  }
}
\keyword{classes}
