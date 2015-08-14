\name{ExpectCovar-class}
\docType{class}
\alias{ExpectCovar-class}
\alias{initialize,ExpectCovar-method}
\title{Class \code{"ExpectCovar"}}
\description{
  Objects of class \code{"ExpectCovar"} represent the expectation and covariance
  of the linear statistic.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{     new("ExpectCovar", \dots).}
}
\section{Slots}{
  \describe{
    \item{\code{expectation}:}{
      Object of class \code{"numeric"}.  The expectation of the linear
      statistic.
    }
    \item{\code{covariance}:}{
      Object of class \code{"matrix"}.  The covariance matrix of the linear
      statistic.
    }
    \item{\code{dimension}:}{
      Object of class \code{"integer"}.  The dimension of the linear statistic.
    }
  }
}
\section{Known Subclasses}{
  Class \code{"\linkS4class{ExpectCovarInfluence}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{
      \code{signature(.Object = "ExpectCovar")}: See the documentation for
      \code{\link[methods:new]{initialize}} (in package \pkg{methods}) for
      details.
    }
  }
}
\section{Note}{
  Normally, this class is used internally only.
}
\keyword{classes}
