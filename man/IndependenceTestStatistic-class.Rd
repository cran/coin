\name{IndependenceTestStatistic-class}
\docType{class}
\alias{IndependenceTestStatistic-class}

\title{Class "IndependenceTestStatistic"}
\description{ Representing test statistics}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IndependenceTestStatistic", itp, varonly)}.
}
\section{Slots}{
	 \describe{
    \item{\code{teststatistic}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{standardizedlinearstatistic}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{linearstatistic}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{expectation}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{covariance}:}{Object of class \code{"VarCovar"} ~~ }
    \item{\code{xtrans}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{ytrans}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{xtrafo}:}{Object of class \code{"function"} ~~ }
    \item{\code{ytrafo}:}{Object of class \code{"function"} ~~ }
    \item{\code{x}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{y}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{weights}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{block}:}{Object of class \code{"factor"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{IndependenceLinearStatistic}"}, directly.
Class \code{"\linkS4class{IndependenceTestProblem}"}, by class "IndependenceLinearStatistic", distance 2.
Class \code{"\linkS4class{IndependenceProblem}"}, by class "IndependenceLinearStatistic", distance 3.
}
\keyword{classes}
