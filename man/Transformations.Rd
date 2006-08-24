\name{Transformations}
\alias{trafo}
\alias{id_trafo}
\alias{ansari_trafo}
\alias{fligner_trafo}
\alias{normal_trafo}
\alias{median_trafo}
\alias{consal_trafo}
\alias{maxstat_trafo}
\alias{fmaxstat_trafo}
\alias{logrank_trafo}
\alias{f_trafo}
\title{ Functions for Data Transformations }
\description{
  Rank-transformations for numerical data or dummy codings of factors.
}
\usage{
trafo(data, numeric_trafo = id_trafo, 
    factor_trafo = f_trafo, surv_trafo = logrank_trafo, 
    var_trafo = NULL, block = NULL)
id_trafo(x)
ansari_trafo(x, ties.method = c("mid-ranks", "average-scores"))
fligner_trafo(x, ties.method = c("mid-ranks", "average-scores"))
normal_trafo(x, ties.method = c("mid-ranks", "average-scores"))
median_trafo(x)
consal_trafo(x, ties.method = c("mid-ranks", "average-scores"))
maxstat_trafo(x, minprob = 0.1, maxprob = 0.9)
logrank_trafo(x, ties.method = c("logrank", "HL"))
f_trafo(x)
}
\arguments{
  \item{data}{an object of class \code{data.frame}.}
  \item{numeric_trafo}{a function to by applied to \code{numeric} 
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns.}
  \item{factor_trafo}{a function to by applied to \code{factor}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns (usually a dummy or contrast 
    matrix).}
  \item{surv_trafo}{a function to by applied to 
    elements of class \code{Surv} of \code{data} returning a 
    matrix with \code{nrow(data)} rows and an arbitrary number of columns.}
  \item{var_trafo}{an optional named list of functions to be applied to the
                   corresponding variables in \code{data}.}
  \item{block}{an optional factor those levels are interpreted as blocks.
    \code{trafo} is applied to each level of \code{block} separately.}
  \item{x}{an object of classes \code{numeric}, \code{factor} or
    \code{Surv}.}
  \item{ties.method}{two methods are available to adjust scores for ties.
    Either the score generating function is applied to \code{mid-ranks} 
    or scores, based on random ranks, are averaged \code{average-scores}.
    For ties handling in case of censored data, see \code{\link{surv_test}}.}
  \item{minprob}{a fraction between 0 and 0.5.}
  \item{maxprob}{a fraction between 0.5 and 1.}
}
\details{

  The utility functions documented here are used to define special
  independence tests.

  \code{trafo} applies its arguments to the elements of \code{data}
  according to the classes of the elements.

  \code{id_trafo} is the identity transformation and \code{f_trafo}
  computes dummy matrices for factors.

  \code{ansari_trafo} and \code{fligner_trafo} compute Ansari-Bradley
  or Fligner scores for scale problems.

  \code{normal_trafo}, \code{median_trafo} and \code{consal_trafo}
  implement normal scores, median scores or Conover-Salburg scores
  (see \code{\link{neuropathy}}) for location problems, 
  \code{logrank_trafo} returns logrank scores for censored data.

  A \code{trafo} function with modified default arguments is usually
  feeded into \code{\link{independence_test}} via the \code{xtrafo}
  or \code{ytrafo} arguments.

  Fine tuning (different transformations for different variables) is
  possible by supplying a named list of functions to the \code{var_trafo}
  argument.

}
\value{

  A named matrix with \code{nrow(data)} rows and 
  arbitrary number of columns. User-supplied transformations must
  return a numeric vector or matrix.

}
\examples{

### dummy matrices, 2-sample problem (only one column)
f_trafo(y <- gl(2, 5))

### K-sample problem (K columns)
f_trafo(y <- gl(5, 2))

### normal scores
normal_trafo(x <- rnorm(10))

### and now together
trafo(data.frame(x = x, y = y), numeric_trafo = normal_trafo)

### the same, more flexible when multiple variables are in play
trafo(data.frame(x = x, y = y), var_trafo = list(x = normal_trafo))

### maximally selected statistics
maxstat_trafo(rnorm(10))

### apply transformation blockwise (e.g. for Friedman test)
trafo(data.frame(y = 1:20), numeric_trafo = rank, block = gl(4, 5))

}
\keyword{manip}
