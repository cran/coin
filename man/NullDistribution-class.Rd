\name{NullDistribution-class}
\docType{class}
\alias{NullDistribution-class}
\alias{ApproxNullDistribution-class}
\alias{AsymptNullDistribution-class}
\alias{ExactNullDistribution-class}
\title{Class \code{"NullDistribution"} and its subclasses}
\description{
  Objects of class \code{"NullDistribution"} and its subclasses
  \code{"ApproxNullDistribution"}, \code{"AsymptNullDistribution"} and
  \code{"ExactNullDistribution"} represent the reference distribution.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{     new("NullDistribution", \dots),

     new("ApproxNullDistribution", \dots),

     new("AsymptNullDistribution", \dots)}
  and \preformatted{     new("ExactNullDistribution", \dots).}
}
\section{Slots}{
  For objects of classes \code{"NullDistribution"},
  \code{"ApproxNullDistribution"}, \code{"AsymptNullDistribution"} or
  \code{"ExactNullDistribution"}:
  \describe{
    \item{\code{q}:}{
      Object of class \code{"function"}.  The quantile function of the reference
      distribution.
    }
    \item{\code{d}:}{
      Object of class \code{"function"}.  The density function of the reference
      distribution.
    }
    \item{\code{support}:}{
      Object of class \code{"function"}.  The support of the reference
      distribution.
    }
    \item{\code{parameters}:}{
      Object of class \code{"list"}.  Additional parameters.
    }
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
  Additionally, for objects of classes \code{"ApproxNullDistribution"} or
  \code{"AsymptNullDistribution"}:
  \describe{
    \item{\code{seed}:}{
      Object of class \code{"integer"}.  The random number generator state
      (i.e., the value of \code{.Random.seed}).
    }
  }
}
\section{Extends}{
  For objects of class \code{"NullDistribution"}: \cr
  Class \code{"\linkS4class{PValue}"}, directly.

  For objects of classes \code{"ApproxNullDistribution"},
  \code{"AsymptNullDistribution"} or \code{"ExactNullDistribution"}: \cr
  Class \code{"NullDistribution"}, directly. \cr
  Class \code{"\linkS4class{PValue}"}, by class \code{"NullDistribution"},
  distance 2.
}
\section{Known Subclasses}{
  For objects of class \code{"NullDistribution"}: \cr
  Class \code{"ApproxNullDistribution"}, directly. \cr
  Class \code{"AsymptNullDistribution"}, directly. \cr
  Class \code{"ExactNullDistribution"}, directly.
}
\section{Methods}{
  \describe{
    \item{dperm}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{dperm}} for details.
    }
    \item{midpvalue}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{midpvalue}} for details.
    }
    \item{pperm}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{pperm}} for details.
    }
    \item{pvalue}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{pvalue}} for details.
    }
    \item{pvalue_interval}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{pvalue_interval}} for details.
    }
    \item{qperm}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{qperm}} for details.
    }
    \item{rperm}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{rperm}} for details.
    }
    \item{support}{
      \code{signature(object = "NullDistribution")}: See the documentation for
      \code{\link{support}} for details.
    }
  }
}
\keyword{classes}
