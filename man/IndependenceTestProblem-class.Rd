\name{IndependenceTestProblem-class}
\docType{class}
\alias{IndependenceTestProblem-class}
\alias{initialize,IndependenceTestProblem-method}

\title{Class "IndependenceTestProblem"}
\description{ Representing an independence problem }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IndependenceTestProblem", ip, xtrafo, ytrafo, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{xtrans}:}{transformed x values }
    \item{\code{ytrans}:}{transformed y values }
    \item{\code{xtrafo}:}{function to transform x values with}
    \item{\code{ytrafo}:}{function to transform y values with }
    \item{\code{x}:}{x values, an object of class \code{"data.frame"}}
    \item{\code{y}:}{y values, an object of class \code{"data.frame"}}
    \item{\code{weights}:}{numeric vector of weights}
    \item{\code{block}:}{factor representing blocks}
  }
}
\section{Extends}{
Class \code{"\linkS4class{IndependenceProblem}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "IndependenceTestProblem")}: ... }
	 }
}
\keyword{classes}
