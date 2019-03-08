\name{MaximallySelectedStatisticsTests}
\alias{maxstat_test}
\alias{maxstat_test.formula}
\alias{maxstat_test.table}
\alias{maxstat_test.IndependenceProblem}
\concept{Generalized maximally selected statistics}
\encoding{UTF-8}
\title{Generalized Maximally Selected Statistics}
\description{
  Testing the independence of two sets of variables measured on arbitrary scales
  against cutpoint alternatives.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{maxstat_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{maxstat_test}{table}(object, \dots)
\method{maxstat_test}{IndependenceProblem}(object, teststat = c("maximum", "quadratic"),
             distribution = c("asymptotic", "approximate", "none"),
             minprob = 0.1, maxprob = 1 - minprob, \dots)
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
    \dots, \code{xp} are considered.  Defaults to \code{0.1}.
  }
  \item{maxprob}{
    a numeric, a fraction between 0.5 and 1 specifying that cutpoints only
    smaller than the \code{maxprob} \eqn{\cdot}{*} 100\% quantile of \code{x1},
    \dots, \code{xp} are considered.  Defaults to \code{1 - minprob}.
  }
  \item{\dots}{
    further arguments to be passed to \code{\link{independence_test}}.
  }
}
\details{
  \code{maxstat_test} provides generalized maximally selected statistics.  The
  family of maximally selected statistics encompasses a large collection of
  procedures used for the estimation of simple cutpoint models including, but
  not limited to, maximally selected \eqn{\chi^2}{chi^2} statistics, maximally
  selected Cochran-Armitage statistics, maximally selected rank statistics and
  maximally selected statistics for multiple covariates.  A general description
  of these methods is given by Hothorn and Zeileis (2008).

  The null hypothesis of independence, or conditional independence given
  \code{block}, between \code{y1}, \dots, \code{yq} and \code{x1}, \dots,
  \code{xp} is tested against cutpoint alternatives.  All possible partitions
  into two groups are evaluated for each unordered covariate \code{x1}, \dots,
  \code{xp}, whereas only order-preserving binary partitions are evaluated for
  ordered or numeric covariates.  The cutpoint is then a set of levels defining
  one of the two groups.

  If both response and covariate is univariable, say \code{y1} and \code{x1},
  this procedure is known as maximally selected \eqn{\chi^2}{chi^2} statistics
  (Miller and Siegmund, 1982) when \code{y1} is a binary factor and \code{x1} is
  a numeric variable, and as maximally selected rank statistics when \code{y1}
  is a rank transformed numeric variable and \code{x1} is a numeric variable
  (Lausen and Schumacher, 1992).  Lausen \emph{et al.} (2004) introduced
  maximally selected statistics for a univariable numeric response and multiple
  numeric covariates \code{x1}, \dots, \code{xp}.

  If, say, \code{y1} and/or \code{x1} are ordered factors, the default scores,
  \code{1:nlevels(y1)} and \code{1:nlevels(x1)} respectively, can be altered
  using the \code{scores} argument (see \code{\link{independence_test}}); this
  argument can also be used to coerce nominal factors to class \code{"ordered"}.
  If both, say, \code{y1} and \code{x1} are ordered factors, a linear-by-linear
  association test is computed and the direction of the alternative hypothesis
  can be specified using the \code{alternative} argument.  The particular
  extension to the case of a univariable binary factor response and a
  univariable ordered covariate was given by Betensky and Rabinowitz (1999) and
  is known as maximally selected Cochran-Armitage statistics.

  The conditional null distribution of the test statistic is used to obtain
  \eqn{p}-values and an asymptotic approximation of the exact distribution is
  used by default (\code{distribution = "asymptotic"}).  Alternatively, the
  distribution can be approximated via Monte Carlo resampling by setting
  \code{distribution} to \code{"approximate"}.  See \code{\link{asymptotic}} and
  \code{\link{approximate}} for details.
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
  Betensky, R. A. and Rabinowitz, D.  (1999).  Maximally selected
  \eqn{\chi^2}{chi^2} statistics for \eqn{k \times 2}{k x 2} tables.
  \emph{Biometrics} \bold{55}(1), 317--320.
  \doi{10.1111/j.0006-341X.1999.00317.x}

  Hothorn, T. and Lausen, B.  (2003).  On the exact distribution of maximally
  selected rank statistics.  \emph{Computational Statistics & Data Analysis}
  \bold{43}(2), 121--137.  \doi{10.1016/S0167-9473(02)00225-6}

  Hothorn, T. and Zeileis, A.  (2008).  Generalized maximally selected
  statistics.  \emph{Biometrics} \bold{64}(4), 1263--1269.
  \doi{10.1111/j.1541-0420.2008.00995.x}

  Lausen, B., Hothorn, T., Bretz, F. and Schumacher, M.  (2004).  Assessment of
  optimal selected prognostic factors.  \emph{Biometrical Journal} \bold{46}(3),
  364--374.  \doi{10.1002/bimj.200310030}

  Lausen, B. and Schumacher, M.  (1992).  Maximally selected rank statistics.
  \emph{Biometrics} \bold{48}(1), 73--85.  \doi{10.2307/2532740}

  Miller, R. and Siegmund, D.  (1982).  Maximally selected chi square
  statistics.  \emph{Biometrics} \bold{38}(4), 1011--1016.
  \doi{10.2307/2529881}

  \enc{Müller}{Mueller}, J. and Hothorn, T.  (2004).  Maximally selected
  two-sample statistics as a new tool for the identification and assessment of
  habitat factors with an application to breeding bird communities in oak
  forests.  \emph{European Journal of Forest Research} \bold{123}(3), 219--228.
  \doi{10.1007/s10342-004-0035-5}
}
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
data("sphase", package = "TH.data")
maxstat_test(Surv(RFS, event) ~  SPF, data = sphase,
             ytrafo = function(data)
                 trafo(data, surv_trafo = function(y)
                     logrank_trafo(y, ties.method = "Hothorn-Lausen")))


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
