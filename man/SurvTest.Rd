\name{SurvTest}
\alias{surv_test}
\alias{surv_test.formula}
\alias{surv_test.IndependenceProblem}
\title{ Independent Two- and K-Sample Tests for Censored Data }
\description{
    
}
\usage{
\method{surv_test}{formula}(formula, data, subset = NULL,  weights = NULL, \dots)
\method{surv_test}{IndependenceProblem}(object, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), ...)
}
\arguments{
  \item{formula}{a formula of the form \code{Surv(time, event) ~ x | block} where 
    \code{time} is a positive numeric variable denoting the survival time and
    \code{event} is a logical being \code{TRUE} when the event of interest
    was observed and \code{FALSE} in case of censoring. \code{x} is a factor
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
    asymptotic distribution (\code{asympt})  
    or via Monte-Carlo resampling (\code{approx}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of the survival
  functions in the groups induced by \code{x} is tested. 

}
\value{

  An object inheriting from class \code{\link{IndependenceTest-class}} with
  methods \code{\link{show}}, \code{\link{statistic}}, \code{\link{expectation}},
  \code{\link{covariance}} and \code{\link{pvalue}}. The null distribution
  can be inspected by \code{\link{pperm}}, \code{\link{dperm}},  
  \code{\link{qperm}} and \code{\link{support}} methods.

}
\examples{

data(ocarcinoma, package = "coin")

surv_test(Surv(time, event) ~ stadium, data = ocarcinoma, 
             distribution = "exact")

}
\keyword{htest}
