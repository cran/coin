\name{MaximallySelectedStatisticsTests}
\alias{maxstat_test}
\alias{maxstat_test.formula}
\alias{maxstat_test.table}
\alias{maxstat_test.IndependenceProblem}
\concept{Generalized maximally selected statistics}
\title{Generalized Maximally Selected Statistics}
\description{
  Testing the independence of two sets of variables measured on arbitrary scales
  against cutpoint alternatives.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{maxstat_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{maxstat_test}{table}(object, ...)
\method{maxstat_test}{IndependenceProblem}(object, teststat = c("maximum", "quadratic"),
             distribution = c("asymptotic", "approximate", "none"),
             minprob = 0.1, maxprob = 1 - minprob, ...)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y1 + ... + yq ~ x1 + ... + xp | block} where
    \code{y1}, \ldots, \code{yq} and \code{x1}, \ldots, \code{xp} are measured
    on arbitrary scales (nominal, ordinal or continuous with or without
    censoring) and \code{block} is an optional factor for stratification.
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
    statistic (\code{"maximum"}, default) or a quadratic form
    (\code{"quadratic"}).
  }
  \item{distribution}{
    a character, the conditional null distribution of the test statistic can be
    approximated by its asymptotic distribution (\code{"asymptotic"}, default)
    or via Monte Carlo resampling (\code{"approximate"}).  Alternatively, the
    functions \code{\link{asymptotic}} or \code{\link{approximate}} can be used.
    Computation of the null distribution can be suppressed by specifying
    \code{"none"}.
  }
  \item{minprob}{
    a numeric, a fraction between 0 and 0.5 specifying that cutpoints only
    greater than the \code{minprob} \eqn{\cdot}{*} 100\% quantile of \code{x1},
    \ldots, \code{xp} are considered.  Defaults to \code{0.1}.
  }
  \item{maxprob}{
    a numeric, a fraction between 0.5 and 1 specifying that cutpoints only
    smaller than the \code{maxprob} \eqn{\cdot}{*} 100\% quantile of \code{x1},
    \ldots, \code{xp} are considered.  Defaults to \code{1 - minprob}.
  }
  \item{...}{
    further arguments to be passed to \code{\link{independence_test}()}.
  }
}
\details{
  \code{maxstat_test()} provides generalized maximally selected statistics.  The
  family of maximally selected statistics encompasses a large collection of
  procedures used for the estimation of simple cutpoint models including, but
  not limited to, maximally selected \eqn{\chi^2}{chi^2} statistics, maximally
  selected Cochran-Armitage statistics, maximally selected rank statistics and
  maximally selected statistics for multiple covariates.  A general description
  of these methods is given by \bibcitet{coin::Hothorn+Zeileis:2008}.

  The null hypothesis of independence, or conditional independence given
  \code{block}, between \code{y1}, \ldots, \code{yq} and \code{x1}, \ldots,
  \code{xp} is tested against cutpoint alternatives.  All possible partitions
  into two groups are evaluated for each unordered covariate \code{x1}, \ldots,
  \code{xp}, whereas only order-preserving binary partitions are evaluated for
  ordered or numeric covariates.  The cutpoint is then a set of levels defining
  one of the two groups.

  If both response and covariate is univariable, say \code{y1} and \code{x1},
  this procedure is known as maximally selected \eqn{\chi^2}{chi^2} statistics
  \bibcitep{coin::miller_1982} when \code{y1} is a binary factor and \code{x1} is
  a numeric variable, and as maximally selected rank statistics when \code{y1}
  is a rank transformed numeric variable and \code{x1} is a numeric variable
  \bibcitep{coin::lausen_1992}. \bibcitet{coin::Lausen:2004:BiomJour} introduced
  maximally selected statistics for a univariable numeric response and multiple
  numeric covariates \code{x1}, \ldots, \code{xp}.

  If, say, \code{y1} and/or \code{x1} are ordered factors, the default scores,
  \code{1:nlevels(y1)} and \code{1:nlevels(x1)}, respectively, can be altered
  using the \code{scores} argument (see \code{\link{independence_test}()}); this
  argument can also be used to coerce nominal factors to class \code{"ordered"}.
  If both, say, \code{y1} and \code{x1} are ordered factors, a linear-by-linear
  association test is computed and the direction of the alternative hypothesis
  can be specified using the \code{alternative} argument.  The particular
  extension to the case of a univariable ordered response and a univariable
  numeric covariate was given by \bibcitet{coin::betensky_1999} and
  is known as maximally selected Cochran-Armitage statistics.

  The conditional null distribution of the test statistic is used to obtain
  \eqn{p}-values and an asymptotic approximation of the exact distribution is
  used by default (\code{distribution = "asymptotic"}).  Alternatively, the
  distribution can be approximated via Monte Carlo resampling by setting
  \code{distribution} to \code{"approximate"}.  See \code{\link{asymptotic}()}
  and \code{\link{approximate}()} for details.

  The tree pipit data used in the example section was described by
  \bibcitet{coin::Mueller:2004:EurJForRes}.
}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\note{
  Starting with version 1.1-0, maximum statistics and quadratic forms can no
  longer be specified using \code{teststat = "maxtype"} and
  \code{teststat = "quadtype"}, respectively (as was used in versions prior to
  0.4-5).
}
\references{\bibshow{*}}
\examples{
\dontshow{options(useFancyQuotes = FALSE)}
## Tree pipit data (Mueller and Hothorn, 2004)
## Asymptotic maximally selected statistics
maxstat_test(counts ~ coverstorey, data = treepipit)

## Asymptotic maximally selected statistics
## Note: all covariates simultaneously
mt <- maxstat_test(counts ~ ., data = treepipit)
mt@estimates$estimate


## Malignant arrythmias data (Hothorn and Lausen, 2003, Sec. 7.2)
## Asymptotic maximally selected statistics
maxstat_test(Surv(time, event) ~  EF, data = hohnloser,
             ytrafo = function(data)
                 trafo(data, surv_trafo = function(y)
                     logrank_trafo(y, ties.method = "Hothorn-Lausen")))


## Breast cancer data (Hothorn and Lausen, 2003, Sec. 7.3)
## Asymptotic maximally selected statistics
if (requireNamespace("TH.data")) {
  data("sphase", package = "TH.data")
  maxstat_test(Surv(RFS, event) ~  SPF, data = sphase,
               ytrafo = function(data)
                   trafo(data, surv_trafo = function(y)
                       logrank_trafo(y, ties.method = "Hothorn-Lausen")))
}


## Job satisfaction data (Agresti, 2002, p. 288, Tab. 7.8)
## Asymptotic maximally selected statistics
maxstat_test(jobsatisfaction)

## Asymptotic maximally selected statistics
## Note: 'Job.Satisfaction' and 'Income' as ordinal
maxstat_test(jobsatisfaction,
             scores = list("Job.Satisfaction" = 1:4,
                           "Income" = 1:4))
}
\keyword{htest}
