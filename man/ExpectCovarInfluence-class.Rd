\name{ExpectCovarInfluence-class}
\docType{class}
\alias{ExpectCovarInfluence-class}

\title{Class "ExpectCovarInfluence"}
\description{  
Conditional expectation and covariance of influence functions
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ExpectCovarInfluence", ...)}.
Normally, this class is used internally only.
}
\section{Slots}{
  \describe{
    \item{\code{sumweights}:}{Object of class \code{"numeric"}.}
    \item{\code{expectation}:}{Object of class \code{"numeric"}.}
    \item{\code{covariance}:}{Object of class \code{"matrix"}.}
    \item{\code{dimension}:}{Object of class \code{"integer"}.}
  }
}
\section{Extends}{
Class \code{"ExpectCovar"}, directly.
}
\keyword{classes}
