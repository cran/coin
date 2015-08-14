\name{PValue-class}
\docType{class}
\alias{PValue-class}
\title{Class \code{"PValue"}}
\description{
  Objects of class \code{"PValue"} represent the \eqn{p}-value,
  mid-\eqn{p}-value and \eqn{p}-value interval of the reference distribution.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{     new("PValue", \dots).}
}
\section{Slots}{
  \describe{
    \item{\code{pvalue}:}{
      Object of class \code{"function"}.  The \eqn{p}-value function of the
      reference distribution.
    }
    \item{\code{midpvalue}:}{
      Object of class \code{"function"}.  The mid-\eqn{p}-value function of the
      reference distribution.
    }
    \item{\code{pvalueinterval}:}{
      Object of class \code{"function"}.  The \eqn{p}-value interval function of
      the reference distribution.
    }
    \item{\code{p}:}{
      Object of class \code{"function"}.  The distribution function of the
      reference distribution.
    }
    \item{\code{name}:}{
      Object of class \code{"character"}.  The name of the reference
      distribution.
    }
  }
}
\section{Methods}{
  No methods defined with class \code{"PValue"} in the signature.
}
\keyword{classes}
