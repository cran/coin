\name{CorrelationTests}
\alias{spearman_test}
\alias{spearman_test.formula}
\alias{spearman_test.IndependenceProblem}
\alias{fisyat_test}
\alias{fisyat_test.formula}
\alias{fisyat_test.IndependenceProblem}
\alias{quadrant_test}
\alias{quadrant_test.formula}
\alias{quadrant_test.IndependenceProblem}
\alias{koziol_test}
\alias{koziol_test.formula}
\alias{koziol_test.IndependenceProblem}
\concept{Spearman correlation test}
\concept{Fisher-Yates correlation test}
\concept{Quadrant test}
\concept{Koziol-Nemec test}
\title{Correlation Tests}
\description{
  Testing the independence of two numeric variables.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{spearman_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{spearman_test}{IndependenceProblem}(object, distribution = c("asymptotic", "approximate", "none"), ...)

\method{fisyat_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{fisyat_test}{IndependenceProblem}(object, distribution = c("asymptotic", "approximate", "none"),
            ties.method = c("mid-ranks", "average-scores"), ...)

\method{quadrant_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{quadrant_test}{IndependenceProblem}(object, distribution = c("asymptotic", "approximate", "none"),
              mid.score = c("0", "0.5", "1"), ...)

\method{koziol_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{koziol_test}{IndependenceProblem}(object, distribution = c("asymptotic", "approximate", "none"),
            ties.method = c("mid-ranks", "average-scores"), ...)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y ~ x | block} where \code{y} and \code{x} are
    numeric variables and \code{block} is an optional factor for stratification.
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
    an object inheriting from class \code{"\linkS4class{IndependenceProblem}"}.
  }
  \item{distribution}{
    a character, the conditional null distribution of the test statistic can be
    approximated by its asymptotic distribution (\code{"asymptotic"}, default)
    or via Monte Carlo resampling (\code{"approximate"}).  Alternatively, the
    functions \code{\link{asymptotic}} or \code{\link{approximate}} can be used.
    Computation of the null distribution can be suppressed by specifying
    \code{"none"}.
  }
  \item{ties.method}{
    a character, the method used to handle ties: the score generating function
    either uses mid-ranks (\code{"mid-ranks"}, default) or averages the scores
    of randomly broken ties (\code{"average-scores"}).
  }
  \item{mid.score}{
    a character, the score assigned to observations exactly equal to the median:
    either 0 (\code{"0"}, default), 0.5 (\code{"0.5"}) or 1 (\code{"1"}); see
    \code{\link{median_test}()}.
  }
  \item{...}{
    further arguments to be passed to \code{\link{independence_test}()}.
  }
}
\details{
  \code{spearman_test()}, \code{fisyat_test()}, \code{quadrant_test()} and
  \code{koziol_test()} provide the Spearman correlation test, the Fisher-Yates
  correlation test using van der Waerden scores, the quadrant test and the
  Koziol-Nemec test.  A general description of these methods is given by
  \enc{Hájek}{Hajek}, \enc{Šidák}{Sidak} and Sen (1999, Sec. 4.6).  The
  Koziol-Nemec test was suggested by Koziol and Nemec (1979).  For the
  adjustment of scores for tied values see \enc{Hájek}{Hajek},
  \enc{Šidák}{Sidak} and Sen (1999, pp. 133--135).

  The null hypothesis of independence, or conditional independence given
  \code{block}, between \code{y} and \code{x} is tested.

  The conditional null distribution of the test statistic is used to obtain
  \eqn{p}-values and an asymptotic approximation of the exact distribution is
  used by default (\code{distribution = "asymptotic"}).  Alternatively, the
  distribution can be approximated via Monte Carlo resampling by setting
  \code{distribution} to \code{"approximate"}.  See \code{\link{asymptotic}()}
  and \code{\link{approximate}()} for details.
}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\references{
  \enc{Hájek}{Hajek}, J., \enc{Šidák}{Sidak}, Z. and Sen, P. K.  (1999).
  \emph{Theory of Rank Tests}, Second Edition.  San Diego: Academic Press.

  Koziol, J. A. and Nemec, A. F.  (1979).  On a \enc{Cramér}{Cramer}-von Mises
  type statistic for testing bivariate independence.  \emph{The Canadian Journal
  of Statistics} \bold{7}(1), 43--52.  \doi{10.2307/3315014}
}
\examples{
## Asymptotic Spearman test
spearman_test(CONT ~ INTG, data = USJudgeRatings)

## Asymptotic Fisher-Yates test
fisyat_test(CONT ~ INTG, data = USJudgeRatings)

## Asymptotic quadrant test
quadrant_test(CONT ~ INTG, data = USJudgeRatings)

## Asymptotic Koziol-Nemec test
koziol_test(CONT ~ INTG, data = USJudgeRatings)
}
\keyword{htest}
