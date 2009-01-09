\name{SurvTest}
\alias{surv_test}
\alias{surv_test.formula}
\alias{surv_test.IndependenceProblem}
\title{ Independent Two- and K-Sample Tests for Censored Data }
\description{
    Testing the equality of survival distributions in
    two or more independent groups.
}
\usage{
\method{surv_test}{formula}(formula, data, subset = NULL,  
    weights = NULL, \dots)
\method{surv_test}{IndependenceProblem}(object, 
    ties.method = c("logrank", "HL", "average-scores"), ...)
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
  \item{ties.method}{a character specifying the way ties are handled in the
                     definition of the logrank scores, see below.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of the survival
  functions in the groups induced by \code{x} is tested. 

  The test implemented here is based on the classical logrank test,
  reformulated as a linear rank test. There are several ways of dealing with
  ties. Here, three methods are implemented. The first one (\code{ties.method = "logrank"})
  is described, for example, in Kalbfleisch \& Prentice (2002, page 221f) or 
  in Callaert (2003) and lead to coefficients 
  \deqn{a_i = \delta_i - \sum_{j: X_j \le X_i} \delta_j / (n - |\{k: X_k < X_j\}|)}
  for a linear rank statistic \eqn{T = \sum_{i = 1}^ n a_i U_i} 
  (in two-sample situations where \eqn{U_i = 0} or \eqn{U_i = 1} denotes
  the groups) with survival times $X_i$ and censoring indicator \eqn{\delta_i = 0}
  for censored observations. For further details, see Kalbfleisch \& Prentice (2002).
  The second method is described in Hothorn \& Lausen (2003) where
  the coefficients
  \deqn{a_i = \delta_i - \sum_{j: X_j \le X_i} \delta_j / (n - |\{k: X_k \le X_j\}| + 1)}
  are suggested. Finally, average scores (as for example used in StatXact) 
  are offered as well.

  Note, however, that the test statistics will differ from the results
  of \code{\link[survival]{survdiff}} since the conditional variance 
  is not identical to the variance estimate used by the classical logrank
  test.

}
\value{

  An object inheriting from class \code{\link{IndependenceTest-class}} with
  methods \code{\link{show}}, \code{\link{statistic}}, \code{\link{expectation}},
  \code{\link{covariance}} and \code{\link{pvalue}}. The null distribution
  can be inspected by \code{\link{pperm}}, \code{\link{dperm}},  
  \code{\link{qperm}} and \code{\link{support}} methods.

}
\references{

    John D. Kalbfleisch \& Ross L. Prentice (2002),
    \emph{The Statistical Analysis of Failure Time Data} (2nd edition).
    John Wiley \& Sons, Hoboken, New Jersey.

    Herman Callaert (2003), Comparing Statistical Software Packages:
    The Case of the Logrank Test in StatXact. 
    \emph{The American Statistician} \bold{57}, 214--217.

    Torsten Hothorn & Berthold Lausen (2003),
    On the Exact Distribution of Maximally Selected Rank
    Statistics. \emph{Computational Statistics \& Data Analysis}
    \bold{43}, 121--137.
}    
\examples{

  ### asymptotic tests for carcinoma data
  surv_test(Surv(time, event) ~ stadium, data = ocarcinoma)
  survdiff(Surv(time, event) ~ stadium, data = ocarcinoma)

  ### example data given in Callaert (2003)
  exdata <- data.frame(time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
                       event = rep(TRUE, 15),
                       group = factor(c(rep(0, 7), rep(1, 8))))
  ### p = 0.0523
  survdiff(Surv(time, event) ~ group, data = exdata)
  ### p = 0.0505
  surv_test(Surv(time, event) ~ group, data = exdata, 
            distribution = exact())
  ### p = 0.0468
  surv_test(Surv(time, event) ~ group, data = exdata, 
            distribution = exact(), ties = "average")

  ### lung cancer example from StatXact
  `lungcancer` <- structure(list(time = c(257, 476, 355, 1779, 355, 191, 
          563, 242, 285, 16, 16, 16, 257, 16), 
      event = c(0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
      group = structure(c(2L, 2L, 2L, 2L, 2L, 1L, 1L, 
          1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("control", "newdrug"), 
          class = "factor")), 
      .Names = c("time", "event", "group"), 
      row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
                    "10", "11", "12", "13", "14"), 
      class = "data.frame")

  ### StatXact 6 manual, page 414
  logrank_trafo(Surv(lungcancer$time, lungcancer$event), 
                ties = "average")

  ### StatXact 6 manual, page 415
  surv_test(Surv(time, event) ~ group, data = lungcancer, 
            ties = "average", distribution = exact())

}
\keyword{htest}
