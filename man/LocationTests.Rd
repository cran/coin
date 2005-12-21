\name{LocationTests}
\alias{oneway_test}
\alias{oneway_test.formula}
\alias{oneway_test.IndependenceProblem}
\alias{wilcox_test.formula}
\alias{wilcox_test.IndependenceProblem}
\alias{wilcox_test}
\alias{normal_test.formula}
\alias{normal_test.IndependenceProblem}
\alias{normal_test}
\alias{median_test.formula}
\alias{median_test.IndependenceProblem}
\alias{median_test}
\alias{kruskal_test.formula}
\alias{kruskal_test.IndependenceProblem}
\alias{kruskal_test}
\title{ Independent Two- and K-Sample Location Tests }
\description{
    Testing the equality of the distributions of a numeric response in
    two or more independent groups against shift alternatives.
}
\usage{

\method{oneway_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{oneway_test}{IndependenceProblem}(object, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"), ...)

\method{wilcox_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{wilcox_test}{IndependenceProblem}(object, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)

\method{normal_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{normal_test}{IndependenceProblem}(object, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"),
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...)

\method{median_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{median_test}{IndependenceProblem}(object, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)

\method{kruskal_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{kruskal_test}{IndependenceProblem}(object, 
    distribution = c("asymptotic", "approximate"), ...)

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a numeric variable giving the data values and \code{x} a factor
    with two or more levels giving the corresponding groups. \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
      integer valued weights for the observations.}
  \item{object}{an object of class \code{IndependenceProblem}.}
  \item{alternative}{a character, the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter.}
  \item{distribution}{a character, the null distribution of the test statistic
    can be computed \code{exact}ly or can be approximated by its
    asymptotic distribution (\code{asymptotic})
    or via Monte-Carlo resampling (\code{approximate}).
    Alternatively, the functions 
    \code{\link{exact}}, \code{\link{approximate}} or \code{\link{asymptotic}} can be
    used to specify how the exact conditional distribution of the test statistic
    should be calculated or approximated.}
  \item{ties.method}{a character, two methods are available to adjust scores for ties,
      either the score generating function is applied to \code{mid-ranks}
      or the scores computed based on random ranks are averaged for all tied
      values (\code{average-scores}).}
  \item{conf.int}{a logical indicating whether a confidence interval
    for the difference in location should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of \code{y} in
  the groups given by \code{x} is tested. In particular, the methods
  documented here are designed to detect shift alternatives. For a general 
  description of the test procedures documented here we refer to Hollander &
  Wolfe (1999).

  The test procedures apply a rank transformation to the response values 
  \code{y}, except of \code{oneway_test} which computes a test statistic
  using the untransformed response values.

  The asymptotic null distribution is computed by default for all
  procedures. Exact p-values may be computed for the two-sample problems and
  can be approximated via Monte-Carlo resampling 
  for all procedures. Exact p-values
  are computed either by the shift algorithm (Streitberg & R\"ohmel, 1986,
  1987) or by the split-up algorithm (van de Wiel, 2001).

  The linear rank tests for two samples (\code{wilcox_test},
  \code{normal_test} and \code{median_test}) can be used to test the 
  two-sided hypothesis \eqn{H_0: Y_1 - Y_2 = 0}, where \eqn{Y_i} is the median
  of the responses in the ith group. Confidence intervals for the difference 
  in location are available for the rank-based procedures and are computed 
  according to Bauer (1972). In case \code{alternative = "less"}, the
  null hypothesis \eqn{H_0: Y_1 - Y_2 \ge 0} is tested and 
  \code{alternative = "greater"} corresponds to a null hypothesis 
  \eqn{H_0: Y_1 - Y_2 \le 0}.

  In case \code{x} is an ordered factor, \code{kruskal_test} computes the
  linear-by-linear association test for ordered alternatives.

  For the adjustment of scores for tied values see Hajek, Sidak and Sen
  (1999), page 131ff.

}
\value{

  An object inheriting from class \code{\link{IndependenceTest-class}} with
  methods \code{\link{show}}, \code{\link{statistic}}, \code{\link{expectation}},
  \code{\link{covariance}} and \code{\link{pvalue}}. The null distribution
  can be inspected by \code{\link{pperm}}, \code{\link{dperm}},
  \code{\link{qperm}} and \code{\link{support}} methods. Confidence
  intervals can be extracted by \code{confint}.

}
\references{

        Myles Hollander & Douglas A. Wolfe (1999),  
        \emph{Nonparametric Statistical Methods, 2nd Edition}. 
        New York: John Wiley & Sons.

        Bernd Streitberg & Joachim R\"ohmel (1986),
        Exact distributions for permutations and rank tests: 
        An introduction to some recently published algorithms.
        \emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

        Bernd Streitberg & Joachim R\"ohmel (1987),
        Exakte Verteilungen f\"ur Rang- und Randomisierungstests
        im allgemeinen $c$-Stichprobenfall.
        \emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

        Mark A. van de Wiel (2001), The split-up algorithm: a fast 
        symbolic method for computing p-values of rank statistics.
        \emph{Computational Statistics} \bold{16}, 519--538.

        David F. Bauer (1972), Constructing confidence sets using 
        rank statistics. \emph{Journal of the American Statistical Association}
        \bold{67}, 687--690.

        Jaroslav Hajek, Zbynek Sidak & Pranab K. Sen (1999), 
        \emph{Theory of Rank Tests}. San Diego, London: Academic Press.

}
\examples{

### Tritiated Water Diffusion Across Human Chorioamnion
### Hollander & Wolfe (1999), Table 4.1, page 110
water_transfer <- data.frame(
    pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
           1.15, 0.88, 0.90, 0.74, 1.21),
    age = factor(c(rep("At term", 10), rep("12-26 Weeks", 5))))

### Wilcoxon-Mann-Whitney test, cf. Hollander & Wolfe (1999), page 111
### exact p-value and confidence interval for the difference in location
### (At term - 12-26 Weeks)
wt <- wilcox_test(pd ~ age, data = water_transfer, 
            distribution = "exact", conf.int = TRUE)
print(wt)

### extract observed Wilcoxon statistic, i.e, the sum of the
### ranks for age = "12-26 Weeks"
statistic(wt, "linear")

### its expectation
expectation(wt)

### and variance
covariance(wt)

### and, finally, the exact two-sided p-value
pvalue(wt)

### Confidence interval for difference (12-26 Weeks - At term)
wilcox_test(pd ~ age, data = water_transfer, 
    xtrafo = function(data) 
        trafo(data, factor_trafo = function(x) as.numeric(x == levels(x)[2])),
        distribution = "exact", conf.int = TRUE)


### Permutation test, asymptotic p-value
oneway_test(pd ~ age, data = water_transfer)

### approximate p-value (with 99\% confidence interval)
pvalue(oneway_test(pd ~ age, data = water_transfer, 
                 distribution = approximate(B = 9999)))
### exact p-value
pt <- oneway_test(pd ~ age, data = water_transfer, distribution = "exact")
pvalue(pt)

### plot density and distribution of the standardised 
### test statistic
layout(matrix(1:2, nrow = 2))
s <- support(pt)
d <- sapply(s, function(x) dperm(pt, x))
p <- sapply(s, function(x) pperm(pt, x))
plot(s, d, type = "S", xlab = "Teststatistic", ylab = "Density")
plot(s, p, type = "S", xlab = "Teststatistic", ylab = "Cumm. Probability")


### Length of YOY Gizzard Shad from Kokosing Lake, Ohio,
### sampled in Summer 1984, Hollander & Wolfe (1999), Table 6.3, page 200
YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44, 
                             42, 60, 32, 42, 45, 58, 27, 51, 42, 52, 
                             38, 33, 26, 25, 28, 28, 26, 27, 27, 27, 
                             31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                  site = factor(c(rep("I", 10), rep("II", 10),
                                  rep("III", 10), rep("IV", 10))))

### Kruskal-Wallis test, approximate exact p-value
kw <- kruskal_test(length ~ site, data = YOY, 
             distribution = approximate(B = 9999))
kw
pvalue(kw)

### Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
### Hollander & Wolfe (1999), page 244 
### (where Steel-Dwass results are given)
if (require("multcomp")) {

    NDWD <- oneway_test(length ~ site, data = YOY,
        ytrafo = function(data) trafo(data, numeric_trafo = rank),
        xtrafo = function(data) trafo(data, factor_trafo = function(x)
            model.matrix(~x - 1) \%*\% t(contrMat(table(x), "Tukey"))),
        teststat = "maxtype", distribution = approximate(B = 90000))

    ### global p-value
    print(pvalue(NDWD))

    ### sites (I = II) != (III = IV) at alpha = 0.01 (page 244)
    print(pvalue(NDWD, method = "single-step"))
}

}
\keyword{htest}
