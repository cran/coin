\name{SymmetryProblem-class}
\docType{class}
\alias{SymmetryProblem-class}
\alias{initialize,SymmetryProblem-method}

\title{Class "SymmetryProblem" }
\description{  The main class for data handling for symmetry problems. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SymmetryProblem", x, y, block, weights)}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"data.frame"}.}
    \item{\code{y}:}{Object of class \code{"data.frame"}.}
    \item{\code{weights}:}{Object of class \code{"numeric"}, however,
      integer values greater or equal zero are expected (and checked).}
    \item{\code{block}:}{Object of class \code{"factor"}.}
  }
}
\section{Extends}{
Class \code{"IndependenceProblem"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "SymmetryProblem")}: ... }
  }
}
\keyword{classes}
