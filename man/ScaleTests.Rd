\name{ScaleTests}
\alias{taha_test}
\alias{taha_test.formula}
\alias{taha_test.IndependenceProblem}
\alias{klotz_test}
\alias{klotz_test.formula}
\alias{klotz_test.IndependenceProblem}
\alias{mood_test}
\alias{mood_test.formula}
\alias{mood_test.IndependenceProblem}
\alias{ansari_test}
\alias{ansari_test.formula}
\alias{ansari_test.IndependenceProblem}
\alias{fligner_test}
\alias{fligner_test.formula}
\alias{fligner_test.IndependenceProblem}
\alias{conover_test}
\alias{conover_test.formula}
\alias{conover_test.IndependenceProblem}
\concept{Taha test}
\concept{Klotz test}
\concept{Mood test}
\concept{Ansari-Bradley test}
\concept{Fligner-Killeen test}
\concept{Conover-Iman test}
\title{Two- and \eqn{K}-Sample Scale Tests}
\description{
  Testing the equality of the distributions of a numeric response variable in
  two or more independent groups against scale alternatives.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{taha_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{taha_test}{IndependenceProblem}(object, conf.int = FALSE, conf.level = 0.95, ...)

\method{klotz_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{klotz_test}{IndependenceProblem}(object, ties.method = c("mid-ranks", "average-scores"),
           conf.int = FALSE, conf.level = 0.95, ...)

\method{mood_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{mood_test}{IndependenceProblem}(object, ties.method = c("mid-ranks", "average-scores"),
          conf.int = FALSE, conf.level = 0.95, ...)

\method{ansari_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{ansari_test}{IndependenceProblem}(object, ties.method = c("mid-ranks", "average-scores"),
            conf.int = FALSE, conf.level = 0.95, ...)

\method{fligner_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{fligner_test}{IndependenceProblem}(object, ties.method = c("mid-ranks", "average-scores"),
             conf.int = FALSE, conf.level = 0.95, ...)

\method{conover_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{conover_test}{IndependenceProblem}(object, conf.int = FALSE, conf.level = 0.95, ...)

}
\arguments{
  \item{formula}{
    a formula of the form \code{y ~ x | block} where \code{y} is a numeric
    variable, \code{x} is a factor and \code{block} is an optional factor for
    stratification.
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
  \item{conf.int}{
    a logical indicating whether a confidence interval for the ratio of scales
    should be computed.  Defaults to \code{FALSE}.
  }
  \item{conf.level}{
    a numeric, confidence level of the interval.  Defaults to \code{0.95}.
  }
  \item{ties.method}{
    a character, the method used to handle ties: the score generating function
    either uses mid-ranks (\code{"mid-ranks"}, default) or averages the scores
    of randomly broken ties (\code{"average-scores"}).
  }
  \item{...}{
    further arguments to be passed to \code{\link{independence_test}()}.
  }
}
\details{
  \code{taha_test()}, \code{klotz_test()}, \code{mood_test()},
  \code{ansari_test()}, \code{fligner_test()} and \code{conover_test()} provide
  the Taha test, the Klotz test, the Mood test, the Ansari-Bradley test, the
  Fligner-Killeen test and the Conover-Iman test.  A general description of
  these methods is given by Hollander and Wolfe (1999).  For the adjustment of
  scores for tied values see \enc{Hájek}{Hajek}, \enc{Šidák}{Sidak} and Sen
  (1999, pp. 133--135).

  The null hypothesis of equality, or conditional equality given \code{block},
  of the distribution of \code{y} in the groups defined by \code{x} is tested
  against scale alternatives.  In the two-sample case, the two-sided null
  hypothesis is \eqn{H_0\!: V(Y_1) / V(Y_2) = 1}{H_0: V(Y_1) / V(Y_2) = 1},
  where \eqn{V(Y_s)} is the variance of the responses in the \eqn{s}th sample.
  In case \code{alternative = "less"}, the null hypothesis is \eqn{H_0\!: V(Y_1)
    / V(Y_2) \ge 1}{H_0: V(Y_1) / V(Y_2) >= 1}.  When
  \code{alternative = "greater"}, the null hypothesis is \eqn{H_0\!: V(Y_1) /
    V(Y_2) \le 1}{H_0: V(Y_1) / V(Y_2) <= 1}.  Confidence intervals for the
  ratio of scales are available and computed according to Bauer (1972).

  The Fligner-Killeen test uses median centering in each of the samples, as
  suggested by Conover, Johnson and Johnson (1981), whereas the Conover-Iman
  test, following Conover and Iman (1978), uses mean centering in each of the
  samples.

  The conditional null distribution of the test statistic is used to obtain
  \eqn{p}-values and an asymptotic approximation of the exact distribution is
  used by default (\code{distribution = "asymptotic"}).  Alternatively, the
  distribution can be approximated via Monte Carlo resampling or computed
  exactly for univariate two-sample problems by setting \code{distribution} to
  \code{"approximate"} or \code{"exact"}, respectively.  See
  \code{\link{asymptotic}()}, \code{\link{approximate}()} and
  \code{\link{exact}()} for details.
}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
  Confidence intervals can be extracted by \code{\link[stats]{confint}()}.
}
\note{
  In the two-sample case, a \emph{large} value of the Ansari-Bradley
  statistic indicates that sample 1 is \emph{less} variable than sample
  2, whereas a \emph{large} value of the statistics due to Taha, Klotz,
  Mood, Fligner-Killeen, and Conover-Iman indicate that sample 1 is
  \emph{more} variable than sample 2.
}
\references{
  Bauer, D. F.  (1972).  Constructing confidence sets using rank statistics.
  \emph{Journal of the American Statistical Association} \bold{67}(339),
  687--690.  \doi{10.1080/01621459.1972.10481279}

  Conover, W. J. and Iman, R. L.  (1978).  Some exact tables for the squared
  ranks test.  \emph{Communications in Statistics -- Simulation and Computation}
  \bold{7}(5), 491--513.  \doi{10.1080/03610917808812093}

  Conover, W. J., Johnson, M. E. and Johnson, M. M.  (1981).  A comparative
  study of tests for homogeneity of variances, with applications to the outer
  continental shelf bidding data.  \emph{Technometrics} \bold{23}(4), 351--361.
  \doi{10.1080/00401706.1981.10487680}

  \enc{Hájek}{Hajek}, J., \enc{Šidák}{Sidak}, Z. and Sen, P. K.  (1999).
  \emph{Theory of Rank Tests}, Second Edition.  San Diego: Academic Press.

  Hollander, M. and Wolfe, D. A.  (1999).  \emph{Nonparametric Statistical
  Methods}, Second Edition.   York: John Wiley & Sons.
}
\examples{
## Serum Iron Determination Using Hyland Control Sera
## Hollander and Wolfe (1999, p. 147, Tab 5.1)
sid <- data.frame(
    serum = c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
              101, 96, 97, 102, 107, 113, 116, 113, 110, 98,
              107, 108, 106, 98, 105, 103, 110, 105, 104,
              100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99),
    method = gl(2, 20, labels = c("Ramsay", "Jung-Parekh"))
)

## Asymptotic Ansari-Bradley test
ansari_test(serum ~ method, data = sid)

## Exact Ansari-Bradley test
pvalue(ansari_test(serum ~ method, data = sid,
                   distribution = "exact"))


## Platelet Counts of Newborn Infants
## Hollander and Wolfe (1999, p. 171, Tab. 5.4)
platelet <- data.frame(
    counts = c(120, 124, 215, 90, 67, 95, 190, 180, 135, 399,
               12, 20, 112, 32, 60, 40),
    treatment = factor(rep(c("Prednisone", "Control"), c(10, 6)))
)

## Approximative (Monte Carlo) Lepage test
## Hollander and Wolfe (1999, p. 172)
lepage_trafo <- function(y)
    cbind("Location" = rank_trafo(y), "Scale" = ansari_trafo(y))

independence_test(counts ~ treatment, data = platelet,
                  distribution = approximate(nresample = 10000),
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = lepage_trafo),
                  teststat = "quadratic")

## Why was the null hypothesis rejected?
## Note: maximum statistic instead of quadratic form
ltm <- independence_test(counts ~ treatment, data = platelet,
                         distribution = approximate(nresample = 10000),
                         ytrafo = function(data)
                             trafo(data, numeric_trafo = lepage_trafo))

## Step-down adjustment suggests a difference in location
pvalue(ltm, method = "step-down")

## The same results are obtained from the simple Sidak-Holm procedure since the
## correlation between Wilcoxon and Ansari-Bradley test statistics is zero
cov2cor(covariance(ltm))
pvalue(ltm, method = "step-down", distribution = "marginal", type = "Sidak")
}
\keyword{htest}
