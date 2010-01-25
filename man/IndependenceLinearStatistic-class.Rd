\name{IndependenceLinearStatistic-class}
\docType{class}
\alias{IndependenceLinearStatistic-class}

\title{Class "IndependenceLinearStatistic"}
\description{ Representing linear statistics }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IndependenceLinearStatistic", ip, xtrafo, ytrafo, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{linearstatistic}:}{Object of class \code{"numeric"}}
    \item{\code{expectation}:}{Object of class \code{"numeric"}}
    \item{\code{covariance}:}{Object of class \code{"VarCovar"}}
    \item{\code{xtrans}:}{Object of class \code{"matrix"}}
    \item{\code{ytrans}:}{Object of class \code{"matrix"}}
    \item{\code{xtrafo}:}{Object of class \code{"function"}}
    \item{\code{ytrafo}:}{Object of class \code{"function"}}
    \item{\code{x}:}{Object of class \code{"data.frame"}}
    \item{\code{y}:}{Object of class \code{"data.frame"}}
    \item{\code{weights}:}{Object of class \code{"numeric"}}
    \item{\code{block}:}{Object of class \code{"factor"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{IndependenceTestProblem}"}, directly.
Class \code{"\linkS4class{IndependenceProblem}"}, by class "IndependenceTestProblem", distance 2.
}
\section{Methods}{
No methods defined with class "IndependenceLinearStatistic" in the signature.
}
\keyword{classes}
