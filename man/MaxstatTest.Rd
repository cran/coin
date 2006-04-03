\name{MaxstatTest}
\alias{maxstat_test}
\alias{maxstat_test.formula}
\alias{maxstat_test.IndependenceProblem}
\title{ Maximally Selected Statistics }
\description{
    Testing the independence of a set of ordered or numeric covariates and a
    response of arbitrary measurement scale against cutpoint alternatives. 
}
\usage{
\method{maxstat_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{maxstat_test}{IndependenceProblem}(object, 
    distribution = c("asymptotic", "approximate"), 
    teststat = c("max", "quad"),
    minprob = 0.1, maxprob = 0.9, ...)
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x1 + ... + xp | block} where \code{y}
    is a variable measured at arbitrary scale and the covariates \code{x1}
    to \code{xp} are at least of class \code{ordered}; \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}   
  \item{object}{an object inheriting from class \code{IndependenceProblem}.}
  \item{distribution}{a character, the null distribution of the test statistic
    can be approximated by its asymptotic distribution (\code{asymptotic}) 
    or via Monte-Carlo resampling (\code{approximate}).
    Alternatively, the functions 
    \code{\link{approximate}} or \code{\link{asymptotic}} can be
    used to specify how the exact conditional distribution of the test statistic
    should be calculated or approximated.}
  \item{teststat}{a character, the type of test statistic to be applied: a
    maximum type statistic (\code{max}) or a quadratic form
    (\code{quad}).}
  \item{minprob}{a fraction between 0 and 0.5; 
    consider only cutpoints greater than 
    the \code{minprob} * 100 \% quantile of \code{x}.}
  \item{maxprob}{a fraction between 0.5 and 1; 
    consider only cutpoints smaller than 
    the \code{maxprob} * 100 \% quantile of \code{x}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of independence of all covariates to the response
  \code{y} against simple cutpoint alternatives is tested.

}
\value{

  An object inheriting from class \code{\link{IndependenceTest-class}} with
  methods \code{\link{show}}, \code{\link{statistic}}, \code{\link{expectation}},
  \code{\link{covariance}} and \code{\link{pvalue}}. The null distribution
  can be inspected by \code{\link{pperm}}, \code{\link{dperm}},  
  \code{\link{qperm}} and \code{\link{support}} methods.

}
\references{

   Rupert Miller & David Siegmund (1982),
   Maximally Selected Chi Square Statistics. 
   \emph{Biometrics} \bold{38}, 1011--1016.

   Berthold Lausen & Martin Schumacher (1992),
   Maximally Selected Rank Statistics.
   \emph{Biometrics} \bold{48}, 73--85.

   Torsten Hothorn & Berthold Lausen (2003),
   On the Exact Distribution of Maximally Selected Rank
   Statistics. \emph{Computational Statistics \& Data Analysis}
   \bold{43}, 121--137.

   Berthold Lausen, Torsten Hothorn, Frank Bretz &
   Martin Schumacher (2004), Optimally Selected Prognostic Factors.
   \emph{Biometrical Journal} \bold{46}, 364--374.

}
\examples{

data("treepipit", package = "coin")

maxstat_test(counts ~ coverstorey, data = treepipit)

}
\keyword{htest}
