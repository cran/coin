\name{pvalue-methods}
\docType{methods}
\alias{pvalue}
\alias{pvalue-methods}
\alias{pvalue,PValue-method}
\alias{pvalue,NullDistribution-method}
\alias{pvalue,ApproxNullDistribution-method}
\alias{pvalue,IndependenceTest-method}
\alias{pvalue,MaxTypeIndependenceTest-method}
\alias{midpvalue}
\alias{midpvalue-methods}
\alias{midpvalue,NullDistribution-method}
\alias{midpvalue,ApproxNullDistribution-method}
\alias{midpvalue,IndependenceTest-method}
\alias{pvalue_interval}
\alias{pvalue_interval-methods}
\alias{pvalue_interval,NullDistribution-method}
\alias{pvalue_interval,IndependenceTest-method}
\alias{size}
\alias{size-methods}
\alias{size,NullDistribution-method}
\alias{size,IndependenceTest-method}
\encoding{UTF-8}
\title{Computation of the \eqn{p}-Value, Mid-\eqn{p}-Value, \eqn{p}-Value
  Interval and Test Size}
\description{
  Methods for computation of the \eqn{p}-value, mid-\eqn{p}-value, \eqn{p}-value
  interval and test size.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\S4method{pvalue}{PValue}(object, q, \dots)
\S4method{pvalue}{NullDistribution}(object, q, \dots)
\S4method{pvalue}{ApproxNullDistribution}(object, q, \dots)
\S4method{pvalue}{IndependenceTest}(object, \dots)
\S4method{pvalue}{MaxTypeIndependenceTest}(object, method = c("global", "single-step",
                          "step-down", "unadjusted"),
       distribution = c("joint", "marginal"),
       type = c("Bonferroni", "Sidak"), \dots)

\S4method{midpvalue}{NullDistribution}(object, q, \dots)
\S4method{midpvalue}{ApproxNullDistribution}(object, q, \dots)
\S4method{midpvalue}{IndependenceTest}(object, \dots)

\S4method{pvalue_interval}{NullDistribution}(object, q, \dots)
\S4method{pvalue_interval}{IndependenceTest}(object, \dots)

\S4method{size}{NullDistribution}(object, alpha, type = c("p-value", "mid-p-value"), \dots)
\S4method{size}{IndependenceTest}(object, alpha, type = c("p-value", "mid-p-value"), \dots)
}
\arguments{
  \item{object}{
    an object from which the \eqn{p}-value, mid-\eqn{p}-value, \eqn{p}-value
    interval or test size can be computed.
  }
  \item{q}{
    a numeric, the quantile for which the \eqn{p}-value, mid-\eqn{p}-value or
    \eqn{p}-value interval is computed.
  }
  \item{method}{
    a character, the method used for the \eqn{p}-value computation: either
    \code{"global"} (default), \code{"single-step"}, \code{"step-down"} or
    \code{"unadjusted"}.
  }
  \item{distribution}{
    a character, the distribution used for the computation of adjusted
    \eqn{p}-values: either \code{"joint"} (default) or \code{"marginal"}.
  }
  \item{type}{
    \code{pvalue()}: a character, the type of \eqn{p}-value adjustment when the
    marginal distributions are used: either \code{"Bonferroni"} (default) or
    \code{"Sidak"}.
    % NOTE: this line is necessary for correct rendering
    \code{size()}: a character, the type of rejection region used when computing
    the test size: either \code{"p-value"} (default) or \code{"mid-p-value"}.
  }
  \item{alpha}{
    a numeric, the nominal significance level \eqn{\alpha} at which the test size
    is computed.
  }
  \item{\dots}{
    further arguments (currently ignored).
  }
}
\details{
  The methods \code{pvalue}, \code{midpvalue}, \code{pvalue_interval} and
  \code{size} compute the \eqn{p}-value, mid-\eqn{p}-value, \eqn{p}-value
  interval and test size respectively.

  For \code{pvalue}, the global \eqn{p}-value (\code{method = "global"}) is
  returned by default and is given with an associated 99\% confidence interval
  when resampling is used to determine the null distribution (which for maximum
  statistics may be true even in the asymptotic case).

  The familywise error rate (FWER) is always controlled under the global null
  hypothesis, i.e., in the \emph{weak} sense, implying that the smallest
  adjusted \eqn{p}-value is valid without further assumptions.  Control of the
  FWER under any partial configuration of the null hypotheses, i.e., in the
  \emph{strong} sense, as is typically desired for multiple tests and
  comparisons, requires that the \emph{subset pivotality} condition holds
  (Westfall and Young, 1993, pp. 42--43; Bretz, Hothorn and Westfall, 2011,
  pp. 136--137).  In addition, for methods based on the joint distribution of
  the test statistics, failure of the \emph{joint exchangeability} assumption
  (Westfall and Troendle, 2008; Bretz, Hothorn and Westfall, 2011, pp. 129--130)
  may cause excess Type I errors.

  Assuming \emph{subset pivotality}, single-step or \emph{free} step-down
  adjusted \eqn{p}-values using max-\eqn{T} procedures are obtained by setting
  \code{method} to \code{"single-step"} or \code{"step-down"} respectively.  In
  both cases, the \code{distribution} argument specifies whether the adjustment
  is based on the joint distribution (\code{"joint"}) or the marginal
  distributions (\code{"marginal"}) of the test statistics.  For procedures
  based on the marginal distributions, Bonferroni- or \enc{Šidák}{Sidak}-type
  adjustment can be specified through the \code{type} argument by setting it to
  \code{"Bonferroni"} or \code{"Sidak"} respectively.

  The \eqn{p}-value adjustment procedures based on the joint distribution of the
  test statistics fully utilizes distributional characteristics, such as
  discreteness and dependence structure, whereas procedures using the marginal
  distributions only incorporate discreteness.  Hence, the joint
  distribution-based procedures are typically more powerful.  Details regarding
  the single-step and \emph{free} step-down procedures based on the joint
  distribution can be found in Westfall and Young (1993); in particular, this
  implementation uses Equation 2.8 with Algorithm 2.5 and 2.8 respectively.
  Westfall and Wolfinger (1997) provide details of the marginal
  distributions-based single-step and \emph{free} step-down procedures.  The
  generalization of Westfall and Wolfinger (1997) to arbitrary test statistics,
  as implemented here, is given by Westfall and Troendle (2008).

  Unadjusted \eqn{p}-values are obtained using \code{method = "unadjusted"}.

  For \code{midpvalue}, the global mid-\eqn{p}-value is given with an associated
  99\% mid-\eqn{p} confidence interval when resampling is used to determine the
  null distribution.  The two-sided mid-\eqn{p}-value is computed according to
  the minimum likelihood method (Hirji \emph{et al.}, 1991).

  The \eqn{p}-value interval \eqn{(p_0, p_1]} obtained by \code{pvalue_interval}
  was proposed by Berger (2000, 2001), where the upper endpoint \eqn{p_1} is the
  conventional \eqn{p}-value and the mid-point, i.e., \eqn{p_{0.5}}{p_0.5}, is
  the mid-\eqn{p}-value.  The lower endpoint \eqn{p_0} is the smallest
  \eqn{p}-value attainable if no conservatism attributable to the discreteness
  of the null distribution is present.  The length of the \eqn{p}-value interval
  is the null probability of the observed outcome and provides a data-dependent
  measure of conservatism that is completely independent of the nominal
  significance level.

  For \code{size}, the test size, i.e., the actual significance level, at the
  nominal significance level \eqn{\alpha} is computed using either the rejection
  region corresponding to the \eqn{p}-value (\code{type = "p-value"}, default)
  or the mid-\eqn{p}-value (\code{type = "mid-p-value"}).  The test size is, in
  contrast to the \eqn{p}-value interval, a data-independent measure of
  conservatism that depends on the nominal significance level.  A test size
  smaller or larger than the nominal significance level indicates that the test
  procedure is conservative or anti-conservative, respectively, at that
  particular nominal significance level.  However, as pointed out by Berger
  (2001), even when the actual and nominal significance levels are identical,
  conservatism may still affect the \eqn{p}-value.
}
\value{
  The \eqn{p}-value, mid-\eqn{p}-value, \eqn{p}-value interval or test size
  computed from \code{object}.  A numeric vector or matrix.
}
\note{
  The mid-\eqn{p}-value, \eqn{p}-value interval and test size of asymptotic
  permutation distributions or exact permutation distributions obtained by the
  split-up algoritm is reported as \code{NA}.

  In versions of \pkg{coin} prior to 1.1-0, a min-\eqn{P} procedure computing
  \enc{Šidák}{Sidak} single-step adjusted \eqn{p}-values accounting for
  discreteness was available when specifying \code{method = "discrete"}.
  \strong{This is now deprecated and will be removed in a future release} due to
  the introduction of a more general max-\eqn{T} version of the same algorithm.
}
\references{
  Berger, V. W.  (2000).  Pros and cons of permutation tests in clinical trials.
  \emph{Statistics in Medicine} \bold{19}(10), 1319--1328.
  \doi{10.1002/(SICI)1097-0258(20000530)19:10<1319::AID-SIM490>3.0.CO;2-0}

  Berger, V. W.  (2001).  The \eqn{p}-value interval as an inferential tool.
  \emph{The Statistician} \bold{50}(1), 79--85.  \doi{10.1111/1467-9884.00262}

  Bretz, F., Hothorn, T. and Westfall, P.  (2011).  \emph{Multiple Comparisons
    Using R}.  Boca Raton: CRC Press.

  Hirji, K. F., Tan, S.-J. and Elashoff, R. M.  (1991).  A quasi-exact test for
  comparing two binomial proportions.  \emph{Statistics in Medicine}
  \bold{10}(7), 1137--1153.  \doi{10.1002/sim.4780100713}

  Westfall, P. H. and Troendle, J. F.  (2008).  Multiple testing with minimal
  assumptions.  \emph{Biometrical Journal} \bold{50}(5), 745--755.
  \doi{10.1002/bimj.200710456}

  Westfall, P. H. and Wolfinger, R. D.  (1997).  Multiple tests with discrete
  distributions.  \emph{The American Statistician} \bold{51}(1), 3--8.
  \doi{10.1080/00031305.1997.10473577}

  Westfall, P. H. and Young, S. S.  (1993).  \emph{Resampling-Based Multiple
    Testing: Examples and Methods for \eqn{p}-Value Adjustment}.  New York: John
  Wiley & Sons.
}
\examples{
## Two-sample problem
dta <- data.frame(
    y = rnorm(20),
    x = gl(2, 10)
)

## Exact Ansari-Bradley test
(at <- ansari_test(y ~ x, data = dta, distribution = "exact"))
pvalue(at)
midpvalue(at)
pvalue_interval(at)
size(at, alpha = 0.05)
size(at, alpha = 0.05, type = "mid-p-value")


## Bivariate two-sample problem
dta2 <- data.frame(
    y1 = rnorm(20) + rep(0:1, each = 10),
    y2 = rnorm(20),
    x = gl(2, 10)
)

## Approximative (Monte Carlo) bivariate Fisher-Pitman test
(it <- independence_test(y1 + y2 ~ x, data = dta2,
                         distribution = approximate(nresample = 10000)))

## Global p-value
pvalue(it)

## Joint distribution single-step p-values
pvalue(it, method = "single-step")

## Joint distribution step-down p-values
pvalue(it, method = "step-down")

## Sidak step-down p-values
pvalue(it, method = "step-down", distribution = "marginal", type = "Sidak")

## Unadjusted p-values
pvalue(it, method = "unadjusted")


## Length of YOY Gizzard Shad (Hollander and Wolfe, 1999, p. 200, Tab. 6.3)
yoy <- data.frame(
    length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
    site = gl(4, 10, labels = as.roman(1:4))
)

## Approximative (Monte Carlo) Fisher-Pitman test with contrasts
## Note: all pairwise comparisons
(it <- independence_test(length ~ site, data = yoy,
                         distribution = approximate(nresample = 10000),
                         xtrafo = mcp_trafo(site = "Tukey")))

## Joint distribution step-down p-values
pvalue(it, method = "step-down") # subset pivotality is violated
}
\keyword{methods}
\keyword{htest}
