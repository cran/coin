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
    minprob = 0.1, maxprob = 1 - minprob, ...)
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x1 + ... + xp | block} where \code{y}
    and covariates \code{x1} to \code{xp} can be variables measured at arbitrary scales; 
    \code{block} is an optional factor for stratification.}
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

  For an unordered covariate \code{x}, all possible partitions into two
  groups are evaluated. The cutpoint is then a set of levels defining 
  one of the two groups.

}
\value{

  An object inheriting from class \code{\link{IndependenceTest-class}} with
  methods \code{\link{show}}, \code{\link{statistic}}, \code{\link{expectation}},
  \code{\link{covariance}} and \code{\link{pvalue}}. The null distribution
  can be inspected by \code{\link{pperm}}, \code{\link{dperm}},  
  \code{\link{qperm}} and \code{\link{support}} methods.

}
\references{

   Rupert Miller \& David Siegmund (1982).
   Maximally Selected Chi Square Statistics. 
   \emph{Biometrics} \bold{38}, 1011--1016.

   Berthold Lausen \& Martin Schumacher (1992).
   Maximally Selected Rank Statistics.
   \emph{Biometrics} \bold{48}, 73--85.

   Torsten Hothorn \& Berthold Lausen (2003).
   On the Exact Distribution of Maximally Selected Rank
   Statistics. \emph{Computational Statistics \& Data Analysis}
   \bold{43}, 121--137.

   Berthold Lausen, Torsten Hothorn, Frank Bretz \&
   Martin Schumacher (2004). Optimally Selected Prognostic Factors.
   \emph{Biometrical Journal} \bold{46}, 364--374.

   J{\"o}rg M{\"u}ller \& Torsten Hothorn (2004).
   Maximally Selected Two-Sample Statistics as a new Tool for
   the Identification and Assessment of Habitat Factors with
   an Application to Breeding Bird Communities in Oak Forests.
   \emph{European Journal of Forest Research}, \bold{123}, 
   218--228.

}
\examples{

### analysis of the tree pipit data in Mueller and Hothorn (2004)
data("treepipit", package = "coin")
maxstat_test(counts ~ coverstorey, data = treepipit)

### and for all possible covariates (simultaneously)
mt <- maxstat_test(counts ~ ., data = treepipit)
show(mt)$estimate

### reproduce applications in Sections 7.2 and 7.3 
### of Hothorn & Lausen (2003) with limiting distribution

data("hohnloser", package = "coin")
maxstat_test(Surv(time, event) ~  EF, data = hohnloser, 
    ytrafo = function(data) trafo(data, surv_trafo = function(x) 
       logrank_trafo(x, ties = "HL")))

data("sphase", package = "coin")
maxstat_test(Surv(RFS, event) ~  SPF, data = sphase,
    ytrafo = function(data) trafo(data, surv_trafo = function(x)
       logrank_trafo(x, ties = "HL")))

}
\keyword{htest}
