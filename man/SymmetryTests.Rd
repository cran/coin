\name{SymmetryTests}
\alias{sign_test}
\alias{sign_test.formula}
\alias{sign_test.SymmetryProblem}
\alias{wilcoxsign_test}
\alias{wilcoxsign_test.formula}
\alias{wilcoxsign_test.SymmetryProblem}
\alias{friedman_test}
\alias{friedman_test.formula}
\alias{friedman_test.SymmetryProblem}
\alias{quade_test}
\alias{quade_test.formula}
\alias{quade_test.SymmetryProblem}
\concept{Sign test}
\concept{Wilcoxon signed-rank test}
\concept{Friedman test}
\concept{Page test}
\concept{Quade test}
\title{Symmetry Tests}
\description{
  Testing the symmetry of a numeric repeated measurements variable in a complete
  block design.
}
\usage{
\method{sign_test}{formula}(formula, data, subset = NULL, \dots)
\method{sign_test}{SymmetryProblem}(object, \dots)

\method{wilcoxsign_test}{formula}(formula, data, subset = NULL, \dots)
\method{wilcoxsign_test}{SymmetryProblem}(object, zero.method = c("Pratt", "Wilcoxon"), \dots)

\method{friedman_test}{formula}(formula, data, subset = NULL, \dots)
\method{friedman_test}{SymmetryProblem}(object, \dots)

\method{quade_test}{formula}(formula, data, subset = NULL, \dots)
\method{quade_test}{SymmetryProblem}(object, \dots)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y ~ x | block} where \code{y} is a numeric
    variable, \code{x} is a factor with two (\code{sign_test} and
    \code{wilcoxsign_test}) or more levels and \code{block} is an optional
    factor (which is generated automatically if omitted).
  }
  \item{data}{
    an optional data frame containing the variables in the model formula.
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be used.  Defaults
    to \code{NULL}.
  }
  \item{object}{
    an object inheriting from class \code{"\linkS4class{SymmetryProblem}"}.
  }
  \item{zero.method}{
    a character, the method used to handle zeros: either \code{"Pratt"}
    (default) or \code{"Wilcoxon"}; see \sQuote{Details}.
  }
  \item{\dots}{
    further arguments to be passed to \code{\link{symmetry_test}}.
  }
}
\details{
  \code{sign_test}, \code{wilcoxsign_test}, \code{friedman_test} and
  \code{quade_test} provide the sign test, the Wilcoxon signed-rank test, the
  Friedman test, the Page test and the Quade test.  A general description of
  these methods is given by Hollander and Wolfe (1999).

  The null hypothesis of symmetry is tested.  The response variable and the
  measurement conditions are given by \code{y} and \code{x}, respectively, and
  \code{block} is a factor where each level corresponds to exactly one subject
  with repeated measurements.  For \code{sign_test} and \code{wilcoxsign_test},
  formulae of the form \code{y ~ x | block} and \code{y ~ x} are allowed.  The
  latter form is interpreted as \code{y} is the first and \code{x} the second
  measurement on the same subject.

  If \code{x} is an ordered factor, the default scores, \code{1:nlevels(x)}, can
  be altered using the \code{scores} argument (see \code{\link{symmetry_test}});
  this argument can also be used to coerce nominal factors to class
  \code{"ordered"}.  In this case, a linear-by-linear association test is
  computed and the direction of the alternative hypothesis can be specified
  using the \code{alternative} argument.  For the Friedman test, this extension
  was given by Page (1963) and is known as the Page test.

  For \code{wilcoxsign_test}, the default method of handling zeros
  (\code{zero.method = "Pratt"}), due to Pratt (1959), first rank transforms the
  absolute differences (including zeros) and then discards the ranks
  corresponding to the zero-differences.  The proposal by Wilcoxon (1949) first
  discards the zero-differences and then rank transforms the remaining absolute
  differences (\code{zero.method = "Wilcoxon"}).

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
  Starting with \pkg{coin} version 1.0-16, the \code{zero.method} argument
  replaced the (now removed) \code{ties.method} argument.  The current default
  is \code{zero.method = "Pratt"} whereas earlier versions had
  \code{ties.method = "HollanderWolfe"}, which is equivalent to
  \code{zero.method = "Wilcoxon"}.
}
\references{
  Hollander, M. and Wolfe, D. A.  (1999).  \emph{Nonparametric Statistical
    Methods}, Second Edition.  New York: John Wiley & Sons.

  Page, E. B.  (1963).  Ordered hypotheses for multiple treatments: a
  significance test for linear ranks.  \emph{Journal of the American Statistical
    Association} \bold{58}(301), 216--230.

  Pratt, J. W.  (1959).  Remarks on zeros and ties in the Wilcoxon signed rank
  procedures.  \emph{Journal of the American Statistical Association}
  \bold{54}(287), 655--667.

  Quade, D.  (1979).  Using weighted rankings in the analysis of complete blocks
  with additive block effects.  \emph{Journal of the American Statistical
    Association} \bold{74}(367), 680--683.

  Wilcoxon, F.  (1949).  Some rapid approximate statistical procedures.
  \emph{Annals of the New York Academy of Sciences} \bold{52}(6), 808--814.
}
\examples{
## Example data from ?wilcox.test
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

## One-sided exact sign test
(st <- sign_test(y1 ~ y2, distribution = "exact",
                 alternative = "greater"))
midpvalue(st) # mid-p-value

## One-sided exact Wilcoxon signed-rank test
(wt <- wilcoxsign_test(y1 ~ y2, distribution = "exact",
                       alternative = "greater"))
statistic(wt, type = "linear")
midpvalue(wt) # mid-p-value

## Comparison with R's wilcox.test() function
wilcox.test(y1, y2, paired = TRUE, alternative = "greater")


## Data with explicit group and block information
dta <- data.frame(y = c(y1, y2), x = gl(2, length(y1)),
                  block = factor(rep(seq_along(y1), 2)))

## For two samples, the sign test is equivalent to the Friedman test...
sign_test(y ~ x | block, data = dta, distribution = "exact")
friedman_test(y ~ x | block, data = dta, distribution = "exact")

## ...and the signed-rank test is equivalent to the Quade test
wilcoxsign_test(y ~ x | block, data = dta, distribution = "exact")
quade_test(y ~ x | block, data = dta, distribution = "exact")


## Comparison of three methods ("round out", "narrow angle", and "wide angle")
## for rounding first base.
## Hollander and Wolfe (1999, p. 274, Tab. 7.1)
rounding <- data.frame(
    times = c(5.40, 5.50, 5.55,
              5.85, 5.70, 5.75,
              5.20, 5.60, 5.50,
              5.55, 5.50, 5.40,
              5.90, 5.85, 5.70,
              5.45, 5.55, 5.60,
              5.40, 5.40, 5.35,
              5.45, 5.50, 5.35,
              5.25, 5.15, 5.00,
              5.85, 5.80, 5.70,
              5.25, 5.20, 5.10,
              5.65, 5.55, 5.45,
              5.60, 5.35, 5.45,
              5.05, 5.00, 4.95,
              5.50, 5.50, 5.40,
              5.45, 5.55, 5.50,
              5.55, 5.55, 5.35,
              5.45, 5.50, 5.55,
              5.50, 5.45, 5.25,
              5.65, 5.60, 5.40,
              5.70, 5.65, 5.55,
              6.30, 6.30, 6.25),
    methods = factor(rep(1:3, 22),
                     labels = c("Round Out", "Narrow Angle", "Wide Angle")),
    block = gl(22, 3)
)

## Asymptotic Friedman test
friedman_test(times ~ methods | block, data = rounding)

## Parallel coordinates plot
with(rounding, {
    matplot(t(matrix(times, ncol = 3, byrow = TRUE)),
            type = "l", lty = 1, col = 1, ylab = "Time", xlim = c(0.5, 3.5),
            axes = FALSE)
    axis(1, at = 1:3, labels = levels(methods))
    axis(2)
})

## Where do the differences come from?
## Wilcoxon-Nemenyi-McDonald-Thompson test (Hollander and Wolfe, 1999, p. 295)
## Note: all pairwise comparisons
(st <- symmetry_test(times ~ methods | block, data = rounding,
                     ytrafo = function(data)
                         trafo(data, numeric_trafo = rank_trafo,
                               block = rounding$block),
                     xtrafo = mcp_trafo(methods = "Tukey")))

## Simultaneous test of all pairwise comparisons
## Wide Angle vs. Round Out differ (Hollander and Wolfe, 1999, p. 296)
pvalue(st, method = "single-step") # subset pivotality is violated


## Strength Index of Cotton
## Hollander and Wolfe (1999, p. 286, Tab. 7.5)
cotton <- data.frame(
    strength = c(7.46, 7.17, 7.76, 8.14, 7.63,
                 7.68, 7.57, 7.73, 8.15, 8.00,
                 7.21, 7.80, 7.74, 7.87, 7.93),
    potash = ordered(rep(c(144, 108, 72, 54, 36), 3),
                     levels = c(144, 108, 72, 54, 36)),
    block = gl(3, 5)
)

## One-sided asymptotic Page test
friedman_test(strength ~ potash | block, data = cotton, alternative = "greater")

## One-sided approximative (Monte Carlo) Page test
friedman_test(strength ~ potash | block, data = cotton, alternative = "greater",
              distribution = approximate(B = 10000))


## Data from Quade (1979, p. 683)
dta <- data.frame(
    y = c(52, 45, 38,
          63, 79, 50,
          45, 57, 39,
          53, 51, 43,
          47, 50, 56,
          62, 72, 49,
          49, 52, 40),
     x = factor(rep(LETTERS[1:3], 7)),
     b = factor(rep(1:7, each = 3))
)

## Approximative (Monte Carlo) Friedman test
## Quade (1979, p. 683)
friedman_test(y ~ x | b, data = dta,
              distribution = approximate(B = 10000)) # chi^2 = 6.000

## Approximative (Monte Carlo) Quade test
## Quade (1979, p. 683)
(qt <- quade_test(y ~ x | b, data = dta,
                  distribution = approximate(B = 10000))) # W = 8.157

## Comparison with R's quade.test() function
quade.test(y ~ x | b, data = dta)

## quade.test() uses an F-statistic
b <- nlevels(qt@statistic@block)
A <- sum(qt@statistic@y^2)
B <- sum(statistic(qt, "linear")^2) / b
(b - 1) * B / (A - B) # F = 8.3765
}
\keyword{htest}
