\name{IndependenceTestStatistic-class}
\docType{class}
\alias{IndependenceTestStatistic-class}
\alias{MaxTypeIndependenceTestStatistic-class}
\alias{QuadTypeIndependenceTestStatistic-class}
\alias{ScalarIndependenceTestStatistic-class}
\alias{initialize,IndependenceTestStatistic-method}
\alias{initialize,MaxTypeIndependenceTestStatistic-method}
\alias{initialize,QuadTypeIndependenceTestStatistic-method}
\alias{initialize,ScalarIndependenceTestStatistic-method}
\title{Class \code{"IndependenceTestStatistic"} and Its Subclasses}
\description{
  Objects of class \code{"IndependenceTestStatistic"} and its subclasses
  \code{"MaxTypeIndependenceTestStatistic"},
  \code{"QuadTypeIndependenceTestStatistic"} and
  \code{"ScalarIndependenceTestStatistic"} represent the test statistic, the
  linear statistic, and the transformed and original data structures
  corresponding to an independence problem.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Class \code{"IndependenceTestStatistic"} is a \emph{virtual} class, so objects
  cannot be created from it directly.

  Objects can be created by calls of the form \preformatted{     new("MaxTypeIndependenceTestStatistic", object,
         alternative = c("two.sided", "less", "greater"), ...),

     new("QuadTypeIndependenceTestStatistic", object, paired = FALSE, ...)}
  and \preformatted{     new("ScalarIndependenceTestStatistic", object,
         alternative = c("two.sided", "less", "greater"), paired = FALSE, ...)}
  where \code{object} is an object of class
  \code{"IndependenceLinearStatistic"}, \code{alternative} is a character
  specifying the direction of the alternative hypothesis and \code{paired} is a
  logical indicating that paired data have been transformed in such a way that
  the (unstandardized) linear statistic is the sum of the absolute values of the
  positive differences between the paired observations.
}
\section{Slots}{
  For objects of classes \code{"IndependenceTestStatistic"},
  \code{"MaxTypeIndependenceTestStatistic"},
  \code{"QuadTypeIndependenceTestStatistic"} or
  \code{"ScalarIndependenceTestStatistic"}:
  \describe{
    \item{\code{teststatistic}:}{
      Object of class \code{"numeric"}.  The test statistic.
    }
    \item{\code{standardizedlinearstatistic}:}{
      Object of class \code{"numeric"}.  The standardized linear statistic.
    }
    \item{\code{linearstatistic}:}{
      Object of class \code{"matrix"}.  The linear statistic for each block.
    }
    \item{\code{expectation}:}{
      Object of class \code{"matrix"}.  The expectation of the linear statistic
      for each block.
    }
    \item{\code{covariance}:}{
      Object of class \code{"matrix"}.  The lower triangular elements of the
      covariance of the linear statistic for each block.
    }
    \item{\code{xtrans}:}{
      Object of class \code{"matrix"}.  The transformed \code{x}.
    }
    \item{\code{ytrans}:}{
      Object of class \code{"matrix"}.  The transformed \code{y}.
    }
    \item{\code{xtrafo}:}{
      Object of class \code{"function"}.  The regression function for \code{x}.
    }
    \item{\code{ytrafo}:}{
      Object of class \code{"function"}.  The influence function for \code{y}.
    }
    \item{\code{x}:}{
      Object of class \code{"data.frame"}.  The variables \code{x}.
    }
    \item{\code{y}:}{
      Object of class \code{"data.frame"}.  The variables \code{y}.
    }
    \item{\code{block}:}{
      Object of class \code{"factor"}.  The block structure.
    }
    \item{\code{weights}:}{
      Object of class \code{"numeric"}.  The case weights.
    }
  }
  Additionally, for objects of classes \code{"MaxTypeIndependenceTest"} or
  \code{"ScalarIndependenceTest"}:
  \describe{
    \item{\code{alternative}:}{
      Object of class \code{"character"}.  The direction of the alternative
      hypothesis.
    }
  }
  Additionally, for objects of class \code{"QuadTypeIndependenceTest"}:
  \describe{
    \item{\code{covarianceplus}:}{
      Object of class \code{"numeric"}.  The lower triangular elements of the
      Moore-Penrose inverse of the covariance of the linear statistic.
    }
    \item{\code{df}:}{
      Object of class \code{"numeric"}.  The rank of the covariance matrix.
    }
  }
  Additionally, for objects of classes \code{"QuadTypeIndependenceTest"} or
  \code{"ScalarIndependenceTest"}:
  \describe{
    \item{\code{paired}:}{
      Object of class \code{"logical"}.  The indicator for paired test
      statistics.
    }
  }
}
\section{Extends}{
  For objects of class \code{"IndependenceTestStatistic"}: \cr
  Class \code{"\linkS4class{IndependenceLinearStatistic}"}, directly. \cr
  Class \code{"\linkS4class{IndependenceTestProblem}"}, by class
  \code{"\linkS4class{IndependenceLinearStatistic}"}, distance 2. \cr
  Class \code{"\linkS4class{IndependenceProblem}"}, by class
  \code{"\linkS4class{IndependenceLinearStatistic}"}, distance 3.

  For objects of classes \code{"MaxTypeIndependenceTestStatistic"},
  \code{"QuadTypeIndependenceTestStatistic"} or
  \code{"ScalarIndependenceTestStatistic"}: \cr
  Class \code{"\linkS4class{IndependenceTestStatistic}"}, directly. \cr
  Class \code{"\linkS4class{IndependenceLinearStatistic}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 2. \cr
  Class \code{"\linkS4class{IndependenceTestProblem}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 3. \cr
  Class \code{"\linkS4class{IndependenceProblem}"}, by class
  \code{"\linkS4class{IndependenceTestStatistic}"}, distance 4.
}
\section{Known Subclasses}{
  For objects of class \code{"IndependenceTestStatistic"}: \cr
  Class \code{"MaxTypeIndependenceTestStatistic"}, directly. \cr
  Class \code{"QuadTypeIndependenceTestStatistic"}, directly. \cr
  Class \code{"ScalarIndependenceTestStatistic"}, directly.
}
\section{Methods}{
  \describe{
    \item{ApproxNullDistribution}{
      \code{signature(object = "MaxTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link{ApproxNullDistribution}()} for details.
    }
    \item{ApproxNullDistribution}{
      \code{signature(object = "QuadTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link{ApproxNullDistribution}()} for details.
    }
    \item{ApproxNullDistribution}{
      \code{signature(object = "ScalarIndependenceTestStatistic")}: See the
      documentation for \code{\link{ApproxNullDistribution}()} for details.
    }
    \item{AsymptNullDistribution}{
      \code{signature(object = "MaxTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link{AsymptNullDistribution}()} for details.
    }
    \item{AsymptNullDistribution}{
      \code{signature(object = "QuadTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link{AsymptNullDistribution}()} for details.
    }
    \item{AsymptNullDistribution}{
      \code{signature(object = "ScalarIndependenceTestStatistic")}: See the
      documentation for \code{\link{AsymptNullDistribution}()} for details.
    }
    \item{ExactNullDistribution}{
      \code{signature(object = "QuadTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link{ExactNullDistribution}()} for details.
    }
    \item{ExactNullDistribution}{
      \code{signature(object = "ScalarIndependenceTestStatistic")}: See the
      documentation for \code{\link{ExactNullDistribution}()} for details.
    }
    \item{covariance}{
      \code{signature(object = "QuadTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link{covariance}()} for details.
    }
    \item{initialize}{
      \code{signature(.Object = "IndependenceTestStatistic")}: See the
      documentation for \code{\link[methods:new]{initialize}()} (in package
      \pkg{methods}) for details.
    }
    \item{initialize}{
      \code{signature(.Object = "MaxTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link[methods:new]{initialize}()} (in package
      \pkg{methods}) for details.
    }
    \item{initialize}{
      \code{signature(.Object = "QuadTypeIndependenceTestStatistic")}: See the
      documentation for \code{\link[methods:new]{initialize}()} (in package
      \pkg{methods}) for details.
    }
    \item{initialize}{
      \code{signature(.Object = "ScalarIndependenceTestStatistic")}: See the
      documentation for \code{\link[methods:new]{initialize}()} (in package
      \pkg{methods}) for details.
    }
    \item{statistic}{
      \code{signature(object = "IndependenceTestStatistic")}: See the
      documentation for \code{\link{statistic}()} for details.
    }
  }
}
\keyword{classes}
