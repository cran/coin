\name{QuadTypeIndependenceTestStatistic-class}
\docType{class}
\alias{QuadTypeIndependenceTestStatistic-class}
\alias{initialize,QuadTypeIndependenceTestStatistic-method}

\title{Class "QuadTypeIndependenceTestStatistic" ~~~ }
\description{ Represent quadratic forms }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("QuadTypeIndependenceTestStatistic", its, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{covarianceplus}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{df}:}{Object of class \code{"numeric"} ~~ }
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
    \item{initialize}{\code{signature(.Object = "QuadTypeIndependenceTestStatistic")}: ... }
	 }
}
\keyword{classes}
