\name{SymmetryTest}
\alias{symmetry_test}
\alias{symmetry_test.formula}
\alias{symmetry_test.table}
\alias{symmetry_test.SymmetryProblem}
\title{General Symmetry Test}
\description{
  Testing the symmetry of set of repeated measurements variables measured on
  arbitrary scales in a complete block design.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{symmetry_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{symmetry_test}{table}(object, \dots)
\method{symmetry_test}{SymmetryProblem}(object, teststat = c("maximum", "quadratic", "scalar"),
              distribution = c("asymptotic", "approximate",
                               "exact", "none"),
              alternative = c("two.sided", "less", "greater"),
              xtrafo = trafo, ytrafo = trafo, scores = NULL,
              check = NULL, paired = FALSE, \dots)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y1 + ... + yq ~ x | block} where \code{y1},
    \dots, \code{yq} are measured on arbitrary scales (nominal, ordinal or
    continuous with or without censoring), \code{x} is a factor and \code{block}
    is an optional factor (which is generated automatically if omitted).
  }
  \item{data}{
    an optional data frame containing the variables in the model formula.
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be used.  Defaults
    to \code{NULL}.
  }
  \item{weights}{
    an optional formula of the form \code{~ w} defining integer valued case
    weights for each observation.  Defaults to \code{NULL}, implying equal
    weight for all observations.
  }
  \item{object}{
    an object inheriting from classes \code{"table"} (with identical
    \code{dimnames} components) or \code{"\linkS4class{SymmetryProblem}"}.
  }
  \item{teststat}{
    a character, the type of test statistic to be applied: either a maximum
    statistic (\code{"maximum"}, default), a quadratic form (\code{"quadratic"})
    or a standardized scalar test statistic (\code{"scalar"}).
  }
  \item{distribution}{
    a character, the conditional null distribution of the test statistic can be
    approximated by its asymptotic distribution (\code{"asymptotic"}, default)
    or via Monte Carlo resampling (\code{"approximate"}).  Alternatively, the
    functions \code{\link{asymptotic}} or \code{\link{approximate}} can be used.
    For univariate two-sample problems, \code{"exact"} or use of the function
    \code{\link{exact}} computes the exact distribution.  Computation of the
    null distribution can be suppressed by specifying \code{"none"}.  It is also
    possible to specify a function with one argument (an object inheriting from
    \code{"\linkS4class{IndependenceTestStatistic}"}) that returns an object of
    class \code{"\linkS4class{NullDistribution}"}.
  }
  \item{alternative}{
    a character, the alternative hypothesis: either \code{"two.sided"}
    (default), \code{"greater"} or \code{"less"}.
  }
  \item{xtrafo}{
    a function of transformations to be applied to the factor \code{x} supplied
    in \code{formula}; see \sQuote{Details}.  Defaults to \code{\link{trafo}}.
  }
  \item{ytrafo}{
    a function of transformations to be applied to the variables \code{y1},
    \dots, \code{yq} supplied in \code{formula}; see \sQuote{Details}.  Defaults
    to \code{\link{trafo}}.
  }
  \item{scores}{
    a named list of scores to be attached to ordered factors; see
    \sQuote{Details}.  Defaults to \code{NULL}, implying equally spaced scores.
  }
  \item{check}{
    a function to be applied to objects of class
    \code{"\linkS4class{IndependenceTest}"} in order to check for specific
    properties of the data.  Defaults to \code{NULL}.
  }
  \item{paired}{
    a logical, indicating that paired data have been transformed in such a way
    that the (unstandardized) linear statistic is the sum of the absolute values
    of the positive differences between the paired observations.  Defaults to
    \code{FALSE}.
  }
  \item{\dots}{
    further arguments to be passed to or from other methods (currently ignored).
  }
}
\details{
  \code{symmetry_test} provides a general symmetry test for a set of variables
  measured on arbitrary scales.  This function is based on the general framework
  for conditional inference procedures proposed by Strasser and Weber (1999).
  The salient parts of the Strasser-Weber framework are elucidated by Hothorn
  \emph{et al.} (2006) and a thorough description of the software implementation
  is given by Hothorn \emph{et al.} (2008).

  The null hypothesis of symmetry is tested.  The response variables and the
  measurement conditions are given by \code{y1}, \dots, \code{yq} and \code{x},
  respectively, and \code{block} is a factor where each level corresponds to
  exactly one subject with repeated measurements.

  A vector of case weights, e.g., observation counts, can be supplied through
  the \code{weights} argument and the type of test statistic is specified by
  the \code{teststat} argument.  Influence and regression functions, i.e.,
  transformations of \code{y1}, \dots, \code{yq} and \code{x}, are specified by
  the \code{ytrafo} and \code{xtrafo} arguments respectively; see
  \code{\link{trafo}} for the collection of transformation functions currently
  available.  This allows for implementation of both novel and familiar test
  statistics, e.g., the McNemar test, the Cochran \eqn{Q} test, the Wilcoxon
  signed-rank test and the Friedman test.  Furthermore, multivariate extensions
  such as the multivariate Friedman test (Gerig, 1969; Puri and Sen, 1971) can
  be implemented without much effort (see \sQuote{Examples}).

  If, say, \code{y1} and/or \code{x} are ordered factors, the default scores,
  \code{1:nlevels(y1)} and \code{1:nlevels(x)} respectively, can be altered
  using the \code{scores} argument; this argument can also be used to coerce
  nominal factors to class \code{"ordered"}.  For example, when \code{y1} is an
  ordered factor with four levels and \code{x} is a nominal factor with three
  levels, \code{scores = list(y1 = c(1, 3:5), x = c(1:2, 4))} supplies the
  scores to be used.  For ordered alternatives the scores must be monotonic, but
  non-montonic scores are also allowed for testing against, e.g., umbrella
  alternatives.  The length of the score vector must be equal to the number of
  factor levels.

  The conditional null distribution of the test statistic is used to obtain
  \eqn{p}-values and an asymptotic approximation of the exact distribution is
  used by default (\code{distribution = "asymptotic"}).  Alternatively, the
  distribution can be approximated via Monte Carlo resampling or computed
  exactly for univariate two-sample problems by setting \code{distribution} to
  \code{"approximate"} or \code{"exact"} respectively.  See
  \code{\link{asymptotic}}, \code{\link{approximate}} and \code{\link{exact}}
  for details.
}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\note{
  Starting with \pkg{coin} version 1.1-0, maximum statistics and quadratic forms
  can no longer be specified using \code{teststat = "maxtype"} and
  \code{teststat = "quadtype"} respectively (as was used in versions prior to
  0.4-5).
}
\references{
  Gerig, T.  (1969).  A multivariate extension of Friedman's
  \eqn{\chi^2_r}{chi^2_r} test.  \emph{Journal of the American Statistical
    Association} \bold{64}(328), 1595--1608.

  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2006).  A Lego
  system for conditional inference.  \emph{The American Statistician}
  \bold{60}(3), 257--263.

  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2008).
  Implementing a class of permutation tests: the coin package.  \emph{Journal of
    Statistical Software} \bold{28}(8), 1--23.
  \url{http://www.jstatsoft.org/v28/i08/}

  Puri, M. L. and Sen, P. K.  (1971).  \emph{Nonparametric Methods in
    Multivariate Analysis}.  New York: John Wiley & Sons.

  Strasser, H. and Weber, C.  (1999).  On the asymptotic theory of permutation
  statistics.  \emph{Mathematical Methods of Statistics} \bold{8}(2), 220--250.
}
\examples{
## One-sided exact Fisher-Pitman test for paired observations
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
dta <- data.frame(
    y = c(y1, y2),
    x = gl(2, length(y1)),
    block = factor(rep(seq_along(y1), 2))
)

symmetry_test(y ~ x | block, data = dta,
              distribution = "exact", alternative = "greater")

## Alternatively: transform data and set 'paired = TRUE'
delta <- y1 - y2
y <- as.vector(rbind(abs(delta) * (delta >= 0), abs(delta) * (delta < 0)))
x <- factor(rep(0:1, length(delta)), labels = c("pos", "neg"))
block <- gl(length(delta), 2)

symmetry_test(y ~ x | block,
              distribution = "exact", alternative = "greater",
              paired = TRUE)


### Example data
### Gerig (1969, p. 1597)
gerig <- data.frame(
    y1 = c( 0.547, 1.811, 2.561,
            1.706, 2.509, 1.414,
           -0.288, 2.524, 3.310,
            1.417, 0.703, 0.961,
            0.878, 0.094, 1.682,
           -0.680, 2.077, 3.181,
            0.056, 0.542, 2.983,
            0.711, 0.269, 1.662,
           -1.335, 1.545, 2.920,
            1.635, 0.200, 2.065),
    y2 = c(-0.575, 1.840, 2.399,
            1.252, 1.574, 3.059,
           -0.310, 1.553, 0.560,
            0.932, 1.390, 3.083,
            0.819, 0.045, 3.348,
            0.497, 1.747, 1.355,
           -0.285, 0.760, 2.332,
            0.089, 1.076, 0.960,
           -0.349, 1.471, 4.121,
            0.845, 1.480, 3.391),
    x = factor(rep(1:3, 10)),
    b = factor(rep(1:10, each = 3))
)

### Asymptotic multivariate Friedman test
### Gerig (1969, p. 1599)
symmetry_test(y1 + y2 ~ x | b, data = gerig, teststat = "quadratic",
              ytrafo = function(data)
                  trafo(data, numeric_trafo = rank_trafo,
                        block = gerig$b)) # L_n = 17.238

### Asymptotic multivariate Page test
(st <- symmetry_test(y1 + y2 ~ x | b, data = gerig,
                     ytrafo = function(data)
                         trafo(data, numeric_trafo = rank_trafo,
                               block = gerig$b),
                     scores = list(x = 1:3)))
pvalue(st, method = "step-down")
}
\keyword{htest}
