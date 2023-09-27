\name{Transformations}
\alias{id_trafo}
\alias{rank_trafo}
\alias{normal_trafo}
\alias{median_trafo}
\alias{savage_trafo}
\alias{consal_trafo}
\alias{koziol_trafo}
\alias{klotz_trafo}
\alias{mood_trafo}
\alias{ansari_trafo}
\alias{fligner_trafo}
\alias{logrank_trafo}
\alias{logrank_weight}
\alias{maxstat_trafo}
\alias{fmaxstat_trafo}
\alias{ofmaxstat_trafo}
\alias{f_trafo}
\alias{of_trafo}
\alias{zheng_trafo}
\alias{trafo}
\alias{mcp_trafo}
\title{Functions for Data Transformation}
\description{
  Transformations for factors and numeric variables.
}
\usage{
id_trafo(x)
rank_trafo(x, ties.method = c("mid-ranks", "random"))
normal_trafo(x, ties.method = c("mid-ranks", "average-scores"))
median_trafo(x, mid.score = c("0", "0.5", "1"))
savage_trafo(x, ties.method = c("mid-ranks", "average-scores"))
consal_trafo(x, ties.method = c("mid-ranks", "average-scores"), a = 5)
koziol_trafo(x, ties.method = c("mid-ranks", "average-scores"), j = 1)
klotz_trafo(x, ties.method = c("mid-ranks", "average-scores"))
mood_trafo(x, ties.method = c("mid-ranks", "average-scores"))
ansari_trafo(x, ties.method = c("mid-ranks", "average-scores"))
fligner_trafo(x, ties.method = c("mid-ranks", "average-scores"))
logrank_trafo(x, ties.method = c("mid-ranks", "Hothorn-Lausen",
                                 "average-scores"),
              weight = logrank_weight, ...)
logrank_weight(time, n.risk, n.event,
               type = c("logrank", "Gehan-Breslow", "Tarone-Ware",
                        "Peto-Peto", "Prentice", "Prentice-Marek",
                        "Andersen-Borgan-Gill-Keiding", "Fleming-Harrington",
                        "Gaugler-Kim-Liao", "Self"),
               rho = NULL, gamma = NULL)
f_trafo(x)
of_trafo(x, scores = NULL)
zheng_trafo(x, increment = 0.1)
maxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
fmaxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
ofmaxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
trafo(data, numeric_trafo = id_trafo, factor_trafo = f_trafo,
      ordered_trafo = of_trafo, surv_trafo = logrank_trafo,
      var_trafo = NULL, block = NULL)
mcp_trafo(...)
}
\arguments{
  \item{x}{
    an object of class \code{"numeric"}, \code{"factor"}, \code{"ordered"} or
    \code{"Surv"}.
  }
  \item{ties.method}{
    a character, the method used to handle ties.  The score generating function
    either uses the mid-ranks (\code{"mid-ranks"}, default) or, in the case of
    \code{rank_trafo()}, randomly broken ties (\code{"random"}).  Alternatively,
    the average of the scores resulting from applying the score generating
    function to randomly broken ties are used (\code{"average-scores"}).  See
    \code{\link{logrank_test}()} for a detailed description of the methods used
    in \code{logrank_trafo()}.
  }
  \item{mid.score}{
    a character, the score assigned to observations exactly equal to the median:
    either 0 (\code{"0"}, default), 0.5 (\code{"0.5"}) or 1 (\code{"1"}); see
    \code{\link{median_test}()}.
  }
  \item{a}{
    a numeric vector, the values taken as the constant \eqn{a} in the
    Conover-Salsburg scores.  Defaults to \code{5}.
  }
  \item{j}{
    a numeric, the value taken as the constant \eqn{j} in the Koziol-Nemec
    scores.  Defaults to \code{1}.
  }
  \item{weight}{
    a function where the first three arguments must correspond to \code{time},
    \code{n.risk}, and \code{n.event} given below.  Defaults to
    \code{logrank_weight}.
  }
  \item{time}{
    a numeric vector, the ordered distinct time points.
  }
  \item{n.risk}{
    a numeric vector, the number of subjects at risk at each time point
    specified in \code{time}.
  }
  \item{n.event}{
    a numeric vector, the number of events at each time point specified in
    \code{time}.
  }
  \item{type}{
    a character, one of \code{"logrank"} (default), \code{"Gehan-Breslow"},
    \code{"Tarone-Ware"}, \code{"Peto-Peto"}, \code{"Prentice"},
    \code{"Prentice-Marek"}, \code{"Andersen-Borgan-Gill-Keiding"},
    \code{"Fleming-Harrington"}, \code{"Gaugler-Kim-Liao"} or \code{"Self"}; see
    \code{\link{logrank_test}()}.
  }
  \item{rho}{
    a numeric vector, the \eqn{\rho} constant when \code{type} is
    \code{"Tarone-Ware"}, \code{"Fleming-Harrington"}, \code{"Gaugler-Kim-Liao"}
    or \code{"Self"}; see \code{\link{logrank_test}()}.  Defaults to
    \code{NULL}, implying \code{0.5} for \code{type = "Tarone-Ware"} and
    \code{0} otherwise.
  }
  \item{gamma}{
    a numeric vector, the \eqn{\gamma} constant when \code{type} is
    \code{"Fleming-Harrington"}, \code{"Gaugler-Kim-Liao"} or \code{"Self"}; see
    \code{\link{logrank_test}()}.  Defaults to \code{NULL}, implying \code{0}.
  }
  \item{scores}{
    a numeric vector or list, the scores corresponding to each level of an
    ordered factor.  Defaults to \code{NULL}, implying \code{1:nlevels(x)}.
  }
  \item{increment}{
    a numeric, the score increment between the order-restricted sets of scores.
    A fraction greater than 0, but smaller than or equal to 1.  Defaults to
    \code{0.1}.
  }
  \item{minprob}{
    a numeric, a fraction between 0 and 0.5; see \code{\link{maxstat_test}()}.
    Defaults to \code{0.1}.
  }
  \item{maxprob}{
    a numeric, a fraction between 0.5 and 1; see \code{\link{maxstat_test}()}.
    Defaults to \code{1 - minprob}.
  }
  \item{data}{
    an object of class \code{"data.frame"}.
  }
  \item{numeric_trafo}{
    a function to be applied to elements of class \code{"numeric"} in
    \code{data}, returning a matrix with \code{nrow(data)} rows and an arbitrary
    number of columns.  Defaults to \code{id_trafo}.
  }
  \item{factor_trafo}{
    a function to be applied to elements of class \code{"factor"} in
    \code{data}, returning a matrix with \code{nrow(data)} rows and an arbitrary
    number of columns.  Defaults to \code{f_trafo}.
  }
  \item{ordered_trafo}{
    a function to be applied to elements of class \code{"ordered"} in
    \code{data}, returning a matrix with \code{nrow(data)} rows and an arbitrary
    number of columns.  Defaults to \code{of_trafo}.
  }
  \item{surv_trafo}{
    a function to be applied to elements of class \code{"Surv"} in \code{data},
    returning a matrix with \code{nrow(data)} rows and an arbitrary number of
    columns.  Defaults to \code{logrank_trafo}.
  }
  \item{var_trafo}{
    an optional named list of functions to be applied to the corresponding
    variables in \code{data}.  Defaults to \code{NULL}.
  }
  \item{block}{
    an optional factor whose levels are interpreted as blocks.  \code{trafo} is
    applied to each level of \code{block} separately.  Defaults to \code{NULL}.
  }
  \item{...}{
    \code{logrank_trafo()}: further arguments to be passed to \code{weight}.
    % NOTE: this line is necessary for correct rendering
    \code{mcp_trafo()}: factor name and contrast matrix (as matrix or character)
    in a \samp{tag = value} format for multiple comparisons based on a single
    unordered factor; see \code{\link[multcomp:glht]{mcp}()} in package
    \pkg{multcomp}.
  }
}
\details{
  The utility functions documented here are used to define specialized test
  procedures.

  \code{id_trafo()} is the identity transformation.

  \code{rank_trafo()}, \code{normal_trafo()}, \code{median_trafo()},
  \code{savage_trafo()}, \code{consal_trafo()} and \code{koziol_trafo()} compute
  rank (Wilcoxon) scores, normal (van der Waerden) scores, median (Mood-Brown)
  scores, Savage scores, Conover-Salsburg scores (see \code{\link{neuropathy}})
  and Koziol-Nemec scores, respectively, for location problems.

  \code{klotz_trafo()}, \code{mood_trafo()}, \code{ansari_trafo()} and
  \code{fligner_trafo()} compute Klotz scores, Mood scores, Ansari-Bradley
  scores and Fligner-Killeen scores, respectively, for scale problems.

  \code{logrank_trafo()} computes weighted logrank scores for right-censored
  data, allowing for a user-defined weight function through the \code{weight}
  argument (see \code{\link{GTSG}}).

  \code{f_trafo()} computes dummy matrices for factors and \code{of_trafo()}
  assigns scores to ordered factors.  For ordered factors with two levels, the
  scores are normalized to the \eqn{[0, 1]} range.  \code{zheng_trafo()}
  computes a finite collection of order-restricted scores for ordered factors
  (see \code{\link{jobsatisfaction}}, \code{\link{malformations}} and
  \code{\link{vision}}).

  \code{maxstat_trafo()}, \code{fmaxstat_trafo()} and \code{ofmaxstat_trafo()}
  compute scores for cutpoint problems (see \code{\link{maxstat_test}()}).

  \code{trafo()} applies its arguments to the elements of \code{data} according
  to the classes of the elements.  A \code{trafo()} function with modified
  default arguments is usually supplied to \code{\link{independence_test}()} via
  the \code{xtrafo} or \code{ytrafo} arguments.  Fine tuning, i.e., different
  transformations for different variables, is possible by supplying a named list
  of functions to the \code{var_trafo} argument.

  \code{mcp_trafo()} computes contrast matrices for factors.
}
\value{
  A numeric vector or matrix with \code{nrow(x)} rows and an arbitrary number of
  columns.  For \code{trafo()}, a named matrix with \code{nrow(data)} rows and an
  arbitrary number of columns.
}
\note{
  Starting with \pkg{coin} version 1.1-0, all transformation functions are now
  passing through missing values (i.e., \code{\link[base]{NA}}s).  Furthermore,
  \code{median_trafo()} and \code{logrank_trafo()} are now \emph{increasing}
  functions (in conformity with most other transformations in this package).
}
\examples{
## Dummy matrix, two-sample problem (only one column)
f_trafo(gl(2, 3))

## Dummy matrix, K-sample problem (K columns)
x <- gl(3, 2)
f_trafo(x)

## Score matrix
ox <- as.ordered(x)
of_trafo(ox)
of_trafo(ox, scores = c(1, 3:4))
of_trafo(ox, scores = list(s1 = 1:3, s2 = c(1, 3:4)))
zheng_trafo(ox, increment = 1/3)

## Normal scores
y <- runif(6)
normal_trafo(y)

## All together now
trafo(data.frame(x = x, ox = ox, y = y), numeric_trafo = normal_trafo)

## The same, but allows for fine-tuning
trafo(data.frame(x = x, ox = ox, y = y), var_trafo = list(y = normal_trafo))

## Transformations for maximally selected statistics
maxstat_trafo(y)
fmaxstat_trafo(x)
ofmaxstat_trafo(ox)

## Apply transformation blockwise (as in the Friedman test)
trafo(data.frame(y = 1:20), numeric_trafo = rank_trafo, block = gl(4, 5))

## Multiple comparisons
dta <- data.frame(x)
mcp_trafo(x = "Tukey")(dta)

## The same, but useful when specific contrasts are desired
K <- rbind("2 - 1" = c(-1,  1, 0),
           "3 - 1" = c(-1,  0, 1),
           "3 - 2" = c( 0, -1, 1))
mcp_trafo(x = K)(dta)
}
\keyword{manip}
