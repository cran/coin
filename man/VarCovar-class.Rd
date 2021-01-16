\name{VarCovar-class}
\docType{class}
\alias{VarCovar-class}
\alias{CovarianceMatrix-class}
\alias{Variance-class}
\alias{initialize,CovarianceMatrix-method}
\alias{initialize,Variance-method}
\title{Class \code{"VarCovar"} and its subclasses}
\description{
  Objects of class \code{"VarCovar"} and its subclasses
  \code{"CovarianceMatrix"} and \code{"Variance"} represent the covariance and
  variance, respectively, of the linear statistic.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Class \code{"VarCovar"} is a \emph{virtual} class defined as the class union
  of \code{"CovarianceMatrix"} and \code{"Variance"}, so objects cannot be
  created from it directly.

  Objects can be created by calls of the form \preformatted{     new("CovarianceMatrix", covariance, \dots)}
  and \preformatted{     new("Variance", variance, \dots)}
  where \code{covariance} is a covariance matrix and \code{variance} is numeric
  vector containing the diagonal elements of the covariance matrix.
}
\section{Slots}{
  For objects of class \code{"CovarianceMatrix"}:
  \describe{
    \item{\code{covariance}:}{
      Object of class \code{"matrix"}.  The covariance matrix.
    }
  }
  For objects of class \code{"Variance"}:
  \describe{
    \item{\code{variance}:}{
      Object of class \code{"numeric"}.  The diagonal elements of the
      covariance matrix.
    }
  }
}
\section{Extends}{
  For objects of classes \code{"CovarianceMatrix"} or \code{"Variance"}: \cr
  Class \code{"VarCovar"}, directly.
}
\section{Known Subclasses}{
  For objects of class \code{"VarCovar"}: \cr
  Class \code{"CovarianceMatrix"}, directly. \cr
  Class \code{"Variance"}, directly.
}
\section{Methods}{
  \describe{
    \item{covariance}{
      \code{signature(object = "CovarianceMatrix")}: See the documentation for
      \code{\link{covariance}} for details.
    }
    \item{initialize}{
      \code{signature(.Object = "CovarianceMatrix")}: See the documentation for
      \code{\link[methods:new]{initialize}} (in package \pkg{methods}) for
      details.
    }
    \item{initialize}{
      \code{signature(.Object = "Variance")}: See the documentation for
      \code{\link[methods:new]{initialize}} (in package \pkg{methods}) for
      details.
    }
    \item{variance}{
      \code{signature(object = "CovarianceMatrix")}: See the documentation for
      \code{\link{variance}} for details.
    }
    \item{variance}{
      \code{signature(object = "Variance")}: See the documentation for
      \code{\link{variance}} for details.
    }
  }
}
\note{
  Starting with \pkg{coin} version 1.4-0, this class is deprecated.  \strong{It
    will be made defunct and removed in a future release.}
}
\keyword{classes}
