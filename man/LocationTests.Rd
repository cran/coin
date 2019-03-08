\name{LocationTests}
\alias{oneway_test}
\alias{oneway_test.formula}
\alias{oneway_test.IndependenceProblem}
\alias{wilcox_test}
\alias{wilcox_test.formula}
\alias{wilcox_test.IndependenceProblem}
\alias{kruskal_test}
\alias{kruskal_test.formula}
\alias{kruskal_test.IndependenceProblem}
\alias{normal_test}
\alias{normal_test.formula}
\alias{normal_test.IndependenceProblem}
\alias{median_test}
\alias{median_test.formula}
\alias{median_test.IndependenceProblem}
\alias{savage_test}
\alias{savage_test.formula}
\alias{savage_test.IndependenceProblem}
\concept{Fisher-Pitman permutation test}
\concept{Wilcoxon-Mann-Whitney test}
\concept{Kruskal-Wallis test}
\concept{van der Waerden test}
\concept{Brown-Mood median test}
\concept{Savage test}
\encoding{UTF-8}
\title{Two- and \eqn{K}-Sample Location Tests}
\description{
  Testing the equality of the distributions of a numeric response variable in
  two or more independent groups against shift alternatives.
}
% NOTE: the markup in the following section is necessary for correct rendering
\usage{
\method{oneway_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{oneway_test}{IndependenceProblem}(object, \dots)

\method{wilcox_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{wilcox_test}{IndependenceProblem}(object, conf.int = FALSE, conf.level = 0.95, \dots)

\method{kruskal_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{kruskal_test}{IndependenceProblem}(object, \dots)

\method{normal_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{normal_test}{IndependenceProblem}(object, ties.method = c("mid-ranks", "average-scores"),
            conf.int = FALSE, conf.level = 0.95, \dots)

\method{median_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{median_test}{IndependenceProblem}(object, mid.score = c("0", "0.5", "1"),
            conf.int = FALSE, conf.level = 0.95, \dots)

\method{savage_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{savage_test}{IndependenceProblem}(object, ties.method = c("mid-ranks", "average-scores"),
            conf.int = FALSE, conf.level = 0.95, \dots)
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
    a logical indicating whether a confidence interval for the difference in
    location should be computed.  Defaults to \code{FALSE}.
  }
  \item{conf.level}{
    a numeric, confidence level of the interval.  Defaults to \code{0.95}.
  }
  \item{ties.method}{
    a character, the method used to handle ties: the score generating function
    either uses mid-ranks (\code{"mid-ranks"}, default) or averages the scores
    of randomly broken ties (\code{"average-scores"}).
  }
  \item{mid.score}{
    a character, the score assigned to observations exactly equal to the median:
    either 0 (\code{"0"}, default), 0.5 (\code{"0.5"}) or 1 (\code{"1"}); see
    \sQuote{Details}.
  }
  \item{\dots}{
    further arguments to be passed to \code{\link{independence_test}}.
  }
}
\details{
  \code{oneway_test}, \code{wilcox_test}, \code{kruskal_test},
  \code{normal_test}, \code{median_test} and \code{savage_test} provide the
  Fisher-Pitman permutation test, the Wilcoxon-Mann-Whitney test, the
  Kruskal-Wallis test, the van der Waerden test, the Brown-Mood median test and
  the Savage test.  A general description of these methods is given by Hollander
  and Wolfe (1999). For the adjustment of scores for tied values see
  \enc{Hájek}{Hajek}, \enc{Šidák}{Sidak} and Sen (1999, pp. 133--135).

  The null hypothesis of equality, or conditional equality given \code{block},
  of the distribution of \code{y} in the groups defined by \code{x} is tested
  against shift alternatives.  In the two-sample case, the two-sided null
  hypothesis is \eqn{H_0\!: \mu = 0}{H_0: mu = 0}, where \eqn{\mu = Y_1 - Y_2}
  and \eqn{Y_s} is the median of the responses in the \eqn{s}th sample.  In case
  \code{alternative = "less"}, the null hypothesis is \eqn{H_0\!: \mu \ge
    0}{H_0: mu >= 0}.  When \code{alternative = "greater"}, the null hypothesis
  is \eqn{H_0\!: \mu \le 0}{H_0: mu <= 0}.  Confidence intervals for the
  difference in location are available (except for \code{oneway_test}) and
  computed according to Bauer (1972).

  If \code{x} is an ordered factor, the default scores, \code{1:nlevels(x)}, can
  be altered using the \code{scores} argument (see
  \code{\link{independence_test}}); this argument can also be used to coerce
  nominal factors to class \code{"ordered"}.  In this case, a linear-by-linear
  association test is computed and the direction of the alternative hypothesis
  can be specified using the \code{alternative} argument.

  The Brown-Mood median test offers a choice of mid-score, i.e., the score
  assigned to observations exactly equal to the median.  In the two-sample case,
  \code{mid-score = "0"} implies that the linear test statistic is simply the
  number of subjects in the second sample with observations greater than the
  median of the pooled sample.  Similarly, the linear test statistic for the
  last alternative, \code{mid-score = "1"}, is the number of subjects in the
  second sample with observations greater than or equal to the median of the
  pooled sample.  If \code{mid-score = "0.5"} is selected, the linear test
  statistic is the mean of the test statistics corresponding to the first and
  last alternatives and has a symmetric distribution, or at least approximately
  so, under the null hypothesis (see \enc{Hájek}{Hajek}, \enc{Šidák}{Sidak}
  and Sen, 1999, pp. 97--98).

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
  Confidence intervals can be extracted by \link[stats]{confint}.
}
\note{
  Starting with version 1.1-0, \code{oneway_test} no longer allows the test
  statistic to be specified; a quadratic form is now used in the \eqn{K}-sample
  case.  Please use \code{\link{independence_test}} if more control is desired.
}
\references{
  Bauer, D. F.  (1972).  Constructing confidence sets using rank statistics.
  \emph{Journal of the American Statistical Association} \bold{67}(339),
  687--690.  \doi{10.1080/01621459.1972.10481279}

  \enc{Hájek}{Hajek}, J., \enc{Šidák}{Sidak}, Z. and Sen, P. K.  (1999).
  \emph{Theory of Rank Tests}, Second Edition.  San Diego: Academic Press.

  Hollander, M. and Wolfe, D. A.  (1999).  \emph{Nonparametric Statistical
    Methods}, Second Edition.  New York: John Wiley & Sons.
}
\examples{\dontshow{options(useFancyQuotes = FALSE)}
## Tritiated Water Diffusion Across Human Chorioamnion
## Hollander and Wolfe (1999, p. 110, Tab. 4.1)
diffusion <- data.frame(
    pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
           1.15, 0.88, 0.90, 0.74, 1.21),
    age = factor(rep(c("At term", "12-26 Weeks"), c(10, 5)))
)

## Exact Wilcoxon-Mann-Whitney test
## Hollander and Wolfe (1999, p. 111)
## (At term - 12-26 Weeks)
(wt <- wilcox_test(pd ~ age, data = diffusion,
                   distribution = "exact", conf.int = TRUE))

## Extract observed Wilcoxon statistic
## Note: this is the sum of the ranks for age = "12-26 Weeks"
statistic(wt, type = "linear")

## Expectation, variance, two-sided pvalue and confidence interval
expectation(wt)
covariance(wt)
pvalue(wt)
confint(wt)

## For two samples, the Kruskal-Wallis test is equivalent to the W-M-W test
kruskal_test(pd ~ age, data = diffusion,
             distribution = "exact")

## Asymptotic Fisher-Pitman test
oneway_test(pd ~ age, data = diffusion)

## Approximative (Monte Carlo) Fisher-Pitman test
pvalue(oneway_test(pd ~ age, data = diffusion,
                   distribution = approximate(nresample = 10000)))

## Exact Fisher-Pitman test
pvalue(ot <- oneway_test(pd ~ age, data = diffusion,
                         distribution = "exact"))

## Plot density and distribution of the standardized test statistic
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:2, nrow = 2))
s <- support(ot)
d <- dperm(ot, s)
p <- pperm(ot, s)
plot(s, d, type = "S", xlab = "Test Statistic", ylab = "Density")
plot(s, p, type = "S", xlab = "Test Statistic", ylab = "Cum. Probability")
par(op) # reset


## Example data
ex <- data.frame(
    y = c(3, 4, 8, 9, 1, 2, 5, 6, 7),
    x = factor(rep(c("no", "yes"), c(4, 5)))
)

## Boxplots
boxplot(y ~ x, data = ex)

## Exact Brown-Mood median test with different mid-scores
(mt1 <- median_test(y ~ x, data = ex, distribution = "exact"))
(mt2 <- median_test(y ~ x, data = ex, distribution = "exact",
                    mid.score = "0.5"))
(mt3 <- median_test(y ~ x, data = ex, distribution = "exact",
                    mid.score = "1")) # sign change!

## Plot density and distribution of the standardized test statistics
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:3, nrow = 3))
s1 <- support(mt1); d1 <- dperm(mt1, s1)
plot(s1, d1, type = "h", main = "Mid-score: 0",
     xlab = "Test Statistic", ylab = "Density")
s2 <- support(mt2); d2 <- dperm(mt2, s2)
plot(s2, d2, type = "h", main = "Mid-score: 0.5",
     xlab = "Test Statistic", ylab = "Density")
s3 <- support(mt3); d3 <- dperm(mt3, s3)
plot(s3, d3, type = "h", main = "Mid-score: 1",
     xlab = "Test Statistic", ylab = "Density")
par(op) # reset


## Length of YOY Gizzard Shad
## Hollander and Wolfe (1999, p. 200, Tab. 6.3)
yoy <- data.frame(
    length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
    site = gl(4, 10, labels = as.roman(1:4))
)

## Approximative (Monte Carlo) Kruskal-Wallis test
kruskal_test(length ~ site, data = yoy,
             distribution = approximate(nresample = 10000))

## Approximative (Monte Carlo) Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
## Hollander and Wolfe (1999, p. 244)
## (where Steel-Dwass results are given)
it <- independence_test(length ~ site, data = yoy,
                        distribution = approximate(nresample = 50000),
                        ytrafo = function(data)
                            trafo(data, numeric_trafo = rank_trafo),
                        xtrafo = mcp_trafo(site = "Tukey"))

## Global p-value
pvalue(it)

## Sites (I = II) != (III = IV) at alpha = 0.01 (p. 244)
pvalue(it, method = "single-step") # subset pivotality is violated
}
\keyword{htest}


%% NOTE:
%% The Jonckhere-Terpstra linear statistic is formulated as the sum of the
%% Mann-Whitney U's for each pairwise comparison, i.e., the outcome is ranked
%% separately for each comparison.  The example below uses the joint ranking and
%% then performs the pairwise comparisons and is *not* equivalent to the J-T
%% test.
%%
%% ## Asymptotic Jonckheere-Terpstra test for ordered groups
%% pieces <- data.frame(
%%     control = c(40, 35, 38, 43, 44, 41),
%%     rough = c(38, 40, 47, 44, 40, 42),
%%     accurate = c(48, 40, 45, 43, 46, 44)
%% )
%% pieces <- stack(pieces)
%% pieces$ind <- ordered(pieces$ind,
%%                       levels = c("control", "rough", "accurate"))
%%
%% ## Look at K: the second line just sums up.
%% ff <- function(x) {
%%     K <- multcomp::contrMat(table(x), "Tukey")[, x]
%%     as.vector(rep(1, nrow(K)) \%*\% K)
%% }
%%
%% independence_test(values ~ ind, data = pieces,
%%                   alternative = "greater",
%%                   ytrafo = function(data)
%%                       trafo(data, numeric_trafo = rank_trafo),
%%                   xtrafo = function(data)
%%                       trafo(data, ordered_trafo = ff))
