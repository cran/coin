\name{ScalarIndependenceTestStatistic-class}
\docType{class}
\alias{ScalarIndependenceTestStatistic-class}
\alias{initialize,ScalarIndependenceTestStatistic-method}

\title{Class "ScalarIndependenceTestStatistic" }
\description{ Represent scalar test statistics }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ScalarIndependenceTestStatistic", its, alternative)}.
}
\section{Slots}{
	 \describe{
    \item{\code{alternative}:}{Object of class \code{"character"} ~~ }
    \item{\code{estimates}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{IndependenceTestStatistic}"}, directly.
Class \code{"\linkS4class{IndependenceLinearStatistic}"}, by class "IndependenceTestStatistic", distance 2.
Class \code{"\linkS4class{IndependenceTestProblem}"}, by class "IndependenceTestStatistic", distance 3.
Class \code{"\linkS4class{IndependenceProblem}"}, by class "IndependenceTestStatistic", distance 4.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ScalarIndependenceTestStatistic")}: ... }
	 }
}
\keyword{classes}
