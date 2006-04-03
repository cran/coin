\name{SymmetryTest}
\alias{symmetry_test}
\alias{symmetry_test.formula}
\alias{symmetry_test.SymmetryProblem}
\alias{symmetry_test.table}
\title{ General Symmetry Test } 
\description{
    Testing the symmetry of a response for repeated measurements in a
    complete block design.
}
\usage{
\method{symmetry_test}{formula}(formula, data, subset = NULL, \dots)
\method{symmetry_test}{SymmetryProblem}(object, 
    teststat = c("max", "quad", "scalar"),
    distribution = c("asymptotic", "approximate", "exact"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL, 
    check = NULL, ...)
\method{symmetry_test}{table}(object, ...)
}
\arguments{
  \item{formula}{a formula of the form 
    \code{y1 + ... + yp ~ x1 + ... + xq | block} where the variables
    on the left and right hand side may be measured on arbitrary scales
    and \code{block} is an optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{object}{an object inheriting from class \code{SymmetryProblem} or an
           object of class \code{table}.}
  \item{teststat}{a character, the type of test statistic to be applied: either a
    standardized scalar test statistic (\code{scalar}), or a
    maximum type statistic (\code{max}) or a quadratic form
    (\code{quad}).}
  \item{alternative}{a character, the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{a character, the null distribution of the test statistic
    can be computed \code{exact}ly or can be approximated by its
    asymptotic distribution (\code{asymptotic})
    or via Monte-Carlo resampling (\code{approximate}).
    Alternatively, the functions 
    \code{\link{exact}}, \code{\link{approximate}} or \code{\link{asymptotic}} can be
    used to specify how the exact conditional distribution of the test statistic
    should be calculated or approximated.}
  \item{xtrafo}{a function of transformations (see \code{\link{trafo}})
    to be applied to the variables on the right hand side of
    \code{formula}, see below.}
  \item{ytrafo}{a function of transformations (see \code{\link{trafo}})
    to be applied to the variables on the left hand side of
    \code{formula}, see below.}
  \item{scores}{a named list of scores to be attached to ordered factors. In
      case a variable is an unordered factor, it is coerced to
      \code{ordered} first.}
  \item{check}{a function to be applied to objects of class
      \code{SymmetryTest} in order to check for specific forms, see below.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  This function handles symmetry problems where each block corresponds to
  exactly one observation with repeated measurements.

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\keyword{htest}
