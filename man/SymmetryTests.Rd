\name{SymmetryTests}
\alias{friedman_test}
\alias{friedman_test.formula}
\alias{friedman_test.SymmetryProblem}
\alias{wilcoxsign_test.formula}
\alias{wilcoxsign_test.IndependenceProblem}
\alias{wilcoxsign_test}
\title{ Symmetry Tests }
\description{
    Testing the symmetry of a response for repeated measurements in a 
    complete block design.
}
\usage{

\method{friedman_test}{formula}(formula, data, subset = NULL, \dots)
\method{friedman_test}{SymmetryProblem}(object, distribution = c("asympt", "approx"), ...) 

\method{wilcoxsign_test}{formula}(formula, data, subset = NULL, \dots)
\method{wilcoxsign_test}{IndependenceProblem}(object, distribution = c("asympt", "approx"), ...) 

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a numeric variable giving the data values and 
    \code{x} a factor with two (\code{wilcoxsign_test}) 
    or more levels giving the corresponding groups. \code{block} is an
    optional factor (which is generated automatically when omitted).}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{object}{an object inheriting from class \code{SymmetryProblem}.}
  \item{distribution}{a character, the null distribution of the test statistic
    can be approximated by its asymptotic distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of symmetry of \code{y} across \code{x} is 
  tested. In case \code{x} is an ordered factor \code{friedman_test}
  performs the Page test, scores can be altered by the \code{scores}
  argument (see \code{\link{symmetry_test}}).

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\references{

        Myles Hollander & Douglas A. Wolfe (1999),  
        \emph{Nonparametric Statistical Methods, 2nd Edition}.
        New York: John Wiley & Sons.
}
\examples{

### Hollander & Wolfe (1999), Table 7.1, page 274
### Comparison of three methods ("round out", "narrow angle", and
###  "wide angle") for rounding first base. 
RoundingTimes <- data.frame(
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
    methods = factor(rep(c("Round Out", "Narrow Angle", "Wide Angle"), 22)),
    block = factor(rep(1:22, rep(3, 22))))

### classical global test
friedman_test(times ~ methods | block, data = RoundingTimes)

### parallel coordinates plot
matplot(t(matrix(RoundingTimes$times, ncol = 3, byrow = TRUE)), 
        type = "l", col = 1, lty = 1, axes = FALSE, ylab = "Time", 
        xlim = c(0.5, 3.5))
axis(1, at = 1:3, labels = levels(RoundingTimes$methods))
axis(2)

### where do the differences come from?
### Wilcoxon-Nemenyi-McDonald-Thompson test
### Hollander & Wolfe (1999), page 295
if (require(multcomp)) {

    ### all pairwise comparisons
    rtt <- symmetry_test(times ~ methods | block, data = RoundingTimes,
         teststat = "maxtype",
         xtrafo = function(data)
             trafo(data, factor_trafo = function(x)
                 model.matrix(~ x - 1) %*% t(contrMat(table(x), "Tukey"))
             ),
         ytrafo = function(data)
             trafo(data, numeric_trafo = rank, block = RoundingTimes$block)
    )

    ### a global test, again
    print(pvalue(rtt))

    ### simultaneous P-values for all pair comparisons
    ### Wide Angle vs. Round Out differ (Hollander and Wolfe, 1999, page 296)
    print(pvalue(rtt, adjusted = TRUE))
}

### Strength Index of Cotton, Hollander & Wolfe (1999), Table 7.5, page 286
sc <- data.frame(block = factor(c(rep(1, 5), rep(2, 5), rep(3, 5))),
                 potash = ordered(rep(c(144, 108, 72, 54, 36), 3)),
                 strength = c(7.46, 7.17, 7.76, 8.14, 7.63,
                              7.68, 7.57, 7.73, 8.15, 8.00,
                              7.21, 7.80, 7.74, 7.87, 7.93))

### Page test for ordered alternatives
ft <- friedman_test(strength ~ potash | block, data = sc)
ft

### one-sided p-value
1 - pnorm(sqrt(statistic(ft)))

### approximate null distribution via Monte-Carlo
pvalue(friedman_test(strength ~ potash | block, data = sc, 
                     distribution = "approx", B = 9999))

}
\keyword{htest}
