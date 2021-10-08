\name{IndependenceTest-class}
\docType{class}
\alias{IndependenceTest-class}
\alias{MaxTypeIndependenceTest-class}
\alias{QuadTypeIndependenceTest-class}
\alias{ScalarIndependenceTest-class}
\alias{ScalarIndependenceTestConfint-class}
\alias{confint,IndependenceTest-method}
\alias{confint,ScalarIndependenceTestConfint-method}
\alias{show,IndependenceTest-method}
\alias{show,MaxTypeIndependenceTest-method}
\alias{show,QuadTypeIndependenceTest-method}
\alias{show,ScalarIndependenceTest-method}
\alias{show,ScalarIndependenceTestConfint-method}
\title{Class \code{"IndependenceTest"} and Its Subclasses}
\description{
  Objects of class \code{"IndependenceTest"} and its subclasses
  \code{"MaxTypeIndependenceTest"}, \code{"QuadTypeIndependenceTest"},
  \code{"ScalarIndependenceTest"} and \code{"ScalarIndependenceTestConfint"}
  represent an independence test including its original and transformed data
  structure, linear statistic, test statistic and reference distribution.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form \preformatted{    new("IndependenceTest", \dots),

    new("MaxTypeIndependenceTest", \dots),

    new("QuadTypeIndependenceTest", \dots),

    new("ScalarIndependenceTest", \dots)}
  and \preformatted{    new("ScalarIndependenceTestConfint", \dots).}
}
\section{Slots}{
  For objects of classes \code{"IndependenceTest"},
  \code{"MaxTypeIndependenceTest"}, \code{"QuadTypeIndependenceTest"},
  \code{"ScalarIndependenceTest"} or \code{"ScalarIndependenceTestConfint"}:
  \describe{
    \item{\code{distribution}:}{
      Object of class \code{"\linkS4class{PValue}"}.  The reference
      distribution.
    }
    \item{\code{statistic}:}{
      Object of class \code{"\linkS4class{IndependenceTestStatistic}"}.  The
      test statistic, the linear statistic, and the transformed and original
      data structures.
    }
    \item{\code{estimates}:}{
      Object of class \code{"list"}.  The estimated parameters.
    }
    \item{\code{method}:}{
      Object of class \code{"character"}.  The test method.
    }
    \item{\code{call}:}{
      Object of class \code{"call"}.  The matched call.
    }
  }
  Additionally, for objects of classes \code{"ScalarIndependenceTest"} or
  \code{"ScalarIndependenceTestConfint"}:
  \describe{
    \item{\code{parameter}:}{
      Object of class \code{"character"}.  The tested parameter.
    }
    \item{\code{nullvalue}:}{
      Object of class \code{"numeric"}.  The hypothesized value of the null
      hypothesis.
    }
  }
  Additionally, for objects of class \code{"ScalarIndependenceTestConfint"}:
  \describe{
    \item{\code{confint}:}{
      Object of class \code{"function"}.  The confidence interval function.
    }
    \item{\code{conf.level}:}{
      Object of class \code{"numeric"}.  The confidence level.
    }
  }
}
\section{Extends}{
  For objects of classes \code{"MaxTypeIndependenceTest"},
  \code{"QuadTypeIndependenceTest"} or \code{"ScalarIndependenceTest"}: \cr
  Class \code{"IndependenceTest"}, directly.

  For objects of class \code{"ScalarIndependenceTestConfint"}: \cr
  Class \code{"ScalarIndependenceTest"}, directly. \cr
  Class \code{"IndependenceTest"}, by class \code{"ScalarIndependenceTest"},
  distance 2.
}
\section{Known Subclasses}{
  For objects of class \code{"IndependenceTest"}: \cr
  Class \code{"MaxTypeIndependenceTest"}, directly. \cr
  Class \code{"QuadTypeIndependenceTest"}, directly. \cr
  Class \code{"ScalarIndependenceTest"}, directly. \cr
  Class \code{"ScalarIndependenceTestConfint"}, by class
  \code{"ScalarIndependenceTest"}, distance 2.

  For objects of class \code{"ScalarIndependenceTest"}: \cr
  Class \code{"ScalarIndependenceTestConfint"}, directly.
}
\section{Methods}{
  \describe{
    \item{confint}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link[stats4]{confint-methods}} (in package \pkg{stats4}) for details.
    }
    \item{confint}{
      \code{signature(object = "ScalarIndependenceTestConfint")}: See the
      documentation for \code{\link[stats4]{confint-methods}} (in package
      \pkg{stats4}) for details.
    }
    \item{covariance}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{covariance}} for details.
    }
    \item{dperm}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{dperm}} for details.
    }
    \item{expectation}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{expectation}} for details.
    }
    \item{midpvalue}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{midpvalue}} for details.
    }
    \item{pperm}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{pperm}} for details.
    }
    \item{pvalue}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{pvalue}} for details.
    }
    \item{pvalue_interval}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{pvalue_interval}} for details.
    }
    \item{qperm}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{qperm}} for details.
    }
    \item{rperm}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{rperm}} for details.
    }
    \item{show}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link[methods]{show}} (in package \pkg{methods}) for details.
    }
    \item{show}{
      \code{signature(object = "MaxTypeIndependenceTest")}: See the
      documentation for \code{\link[methods]{show}} (in package \pkg{methods})
      for details.
    }
    \item{show}{
      \code{signature(object = "QuadTypeIndependenceTest")}: See the
      documentation for \code{\link[methods]{show}} (in package \pkg{methods})
      for details.
    }
    \item{show}{
      \code{signature(object = "ScalarIndependenceTest")}: See the documentation
      for \code{\link[methods]{show}} (in package \pkg{methods}) for details.
    }
    \item{show}{
      \code{signature(object = "ScalarIndependenceTestConfint")}: See the
      documentation for \code{\link[methods]{show}} (in package \pkg{methods})
      for details.
    }
    \item{size}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{size}} for details.
    }
    \item{statistic}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{statistic}} for details.
    }
    \item{support}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{support}} for details.
    }
    \item{variance}{
      \code{signature(object = "IndependenceTest")}: See the documentation for
      \code{\link{variance}} for details.
    }
  }
}
\keyword{classes}
