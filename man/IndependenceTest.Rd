\name{IndependenceTest}
\alias{independence_test}
\alias{independence_test.formula}
\alias{independence_test.table}
\alias{independence_test.IndependenceProblem}
\title{General Independence Test}
\description{
  Testing the independence of two sets of variables measured on arbitrary scales.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{independence_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{independence_test}{table}(object, \dots)
\method{independence_test}{IndependenceProblem}(object, teststat = c("maximum", "quadratic", "scalar"),
                  distribution = c("asymptotic", "approximate",
                                   "exact", "none"),
                  alternative = c("two.sided", "less", "greater"),
                  xtrafo = trafo, ytrafo = trafo, scores = NULL,
                  check = NULL, \dots)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y1 + ... + yq ~ x1 + ... + xp | block} where
    \code{y1}, \dots, \code{yq} and \code{x1}, \dots, \code{xp} are measured on
    arbitrary scales (nominal, ordinal or continuous with or without censoring)
    and \code{block} is an optional factor for stratification.
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
    an object inheriting from classes \code{"table"} or
    \code{"\linkS4class{IndependenceProblem}"}.
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
    a function of transformations to be applied to the variables \code{x1},
    \dots, \code{xp} supplied in \code{formula}; see \sQuote{Details}.  Defaults
    to \code{\link{trafo}}.
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
  \item{\dots}{
    further arguments to be passed to or from other methods (currently ignored).
  }
}
\details{
  \code{independence_test} provides a general independence test for two sets of
  variables measured on arbitrary scales.  This function is based on the general
  framework for conditional inference procedures proposed by Strasser and Weber
  (1999).  The salient parts of the Strasser-Weber framework are elucidated by
  Hothorn \emph{et al.} (2006) and a thorough description of the software
  implementation is given by Hothorn \emph{et al.} (2008).

  The null hypothesis of independence, or conditional independence given
  \code{block}, between \code{y1}, \dots, \code{yq} and \code{x1}, \dots,
  \code{xp} is tested.

  A vector of case weights, e.g., observation counts, can be supplied through
  the \code{weights} argument and the type of test statistic is specified by the
  \code{teststat} argument.  Influence and regression functions, i.e.,
  transformations of \code{y1}, \dots, \code{yq} and \code{x1}, \dots,
  \code{xp}, are specified by the \code{ytrafo} and \code{xtrafo} arguments
  respectively; see \code{\link{trafo}} for the collection of transformation
  functions currently available.  This allows for implementation of both novel
  and familiar test statistics, e.g., the Pearson \eqn{\chi^2} test, the
  generalized Cochran-Mantel-Haenszel test, the Spearman correlation test, the
  Fisher-Pitman permutation test, the Wilcoxon-Mann-Whitney test, the
  Kruskal-Wallis test and the family of weighted logrank tests for censored
  data.  Furthermore, multivariate extensions such as the multivariate
  Kruskal-Wallis test (Puri and Sen, 1966, 1971) can be implemented without much
  effort (see \sQuote{Examples}).

  If, say, \code{y1} and/or \code{x1} are ordered factors, the default scores,
  \code{1:nlevels(y1)} and \code{1:nlevels(x1)} respectively, can be altered
  using the \code{scores} argument; this argument can also be used to coerce
  nominal factors to class \code{"ordered"}.  For example, when \code{y1} is an
  ordered factor with four levels and \code{x1} is a nominal factor with three
  levels, \code{scores = list(y1 = c(1, 3:5), x1 = c(1:2, 4))} supplies the
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
  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2006).  A Lego
  system for conditional inference.  \emph{The American Statistician}
  \bold{60}(3), 257--263.

  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2008).
  Implementing a class of permutation tests: the coin package.  \emph{Journal of
    Statistical Software} \bold{28}(8), 1--23.
  \url{http://www.jstatsoft.org/v28/i08/}

  Johnson, W. D., Mercante, D. E. and May, W. L.  (1993).  A computer package
  for the multivariate nonparametric rank test in completely randomized
  experimental designs.  \emph{Computer Methods and Programs in Biomedicine}
  \bold{40}(3), 217--225.

  Puri, M. L. and Sen, P. K.  (1966).  On a class of multivariate multisample
  rank order tests.  \emph{Sankhya} A \bold{28}(4), 353--376.

  Puri, M. L. and Sen, P. K.  (1971).  \emph{Nonparametric Methods in
    Multivariate Analysis}.  New York: John Wiley & Sons.

  Strasser, H. and Weber, C.  (1999).  On the asymptotic theory of permutation
  statistics.  \emph{Mathematical Methods of Statistics} \bold{8}(2), 220--250.
}
\examples{
## One-sided exact van der Waerden (normal scores) test...
independence_test(asat ~ group, data = asat,
                  ## exact null distribution
                  distribution = "exact",
                  ## one-sided test
                  alternative = "greater",
                  ## apply normal scores to asat$asat
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = normal_trafo),
                  ## indicator matrix of 1st level of asat$group
                  xtrafo = function(data)
                      trafo(data, factor_trafo = function(x)
                          matrix(x == levels(x)[1], ncol = 1)))

## ...or more conveniently
normal_test(asat ~ group, data = asat,
            ## exact null distribution
            distribution = "exact",
            ## one-sided test
            alternative = "greater")


## Receptor binding assay of benzodiazepines
## Johnson, Mercante and May (1993, Tab. 1)
benzos <- data.frame(
      cerebellum = c( 3.41,  3.50,  2.85,  4.43,
                      4.04,  7.40,  5.63, 12.86,
                      6.03,  6.08,  5.75,  8.09,  7.56),
       brainstem = c( 3.46,  2.73,  2.22,  3.16,
                      2.59,  4.18,  3.10,  4.49,
                      6.78,  7.54,  5.29,  4.57,  5.39),
          cortex = c(10.52,  7.52,  4.57,  5.48,
                      7.16, 12.00,  9.36,  9.35,
                     11.54, 11.05,  9.92, 13.59, 13.21),
    hypothalamus = c(19.51, 10.00,  8.27, 10.26,
                     11.43, 19.13, 14.03, 15.59,
                     24.87, 14.16, 22.68, 19.93, 29.32),
        striatum = c( 6.98,  5.07,  3.57,  5.34,
                      4.57,  8.82,  5.76, 11.72,
                      6.98,  7.54,  7.66,  9.69,  8.09),
     hippocampus = c(20.31, 13.20,  8.58, 11.42,
                     13.79, 23.71, 18.35, 38.52,
                     21.56, 18.66, 19.24, 27.39, 26.55),
       treatment = factor(rep(c("Lorazepam", "Alprazolam", "Saline"),
                          c(4, 4, 5)))
)

## Approximative (Monte Carlo) multivariate Kruskal-Wallis test
## Johnson, Mercante and May (1993, Tab. 2)
independence_test(cerebellum + brainstem + cortex +
                  hypothalamus + striatum + hippocampus ~ treatment,
                  data = benzos,
                  teststat = "quadratic",
                  distribution = approximate(B = 10000),
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = rank_trafo)) # Q = 16.129
}
\keyword{htest}
