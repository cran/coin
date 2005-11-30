\name{ContingencyTests}
\alias{chisq_test}
\alias{chisq_test.formula}
\alias{chisq_test.table}
\alias{chisq_test.IndependenceProblem}
\alias{cmh_test.formula}
\alias{cmh_test.table}
\alias{cmh_test.IndependenceProblem}
\alias{cmh_test}
\alias{lbl_test.formula}
\alias{lbl_test.table}
\alias{lbl_test.IndependenceProblem}
\alias{lbl_test}
\title{ Independence in Three-Way Contingency Tables }
\description{
    Testing the independence of two possibly ordered factors, eventually
    stratified by a third factor.
}
\usage{

\method{cmh_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{cmh_test}{table}(object, distribution = c("asymptotic", "approximate"), ...)
\method{cmh_test}{IndependenceProblem}(object, distribution = c("asymptotic", "approximate"), ...)

\method{chisq_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{chisq_test}{table}(object, distribution = c("asymptotic", "approximate"), ...)
\method{chisq_test}{IndependenceProblem}(object, distribution = c("asymptotic", "approximate"), ...)

\method{lbl_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{lbl_test}{table}(object, distribution = c("asymptotic", "approximate"), ...)
\method{lbl_test}{IndependenceProblem}(object, distribution = c("asymptotic", "approximate"), ...)

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    and \code{x} are factors (possibly ordered) and \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}
  \item{object}{an object inheriting from class \code{"IndependenceProblem"} or an
    object of class \code{table}.}
  \item{distribution}{a character, the null distribution of the test statistic
    can be approximated by its asymptotic distribution (\code{"asymptotic"}) 
    or via Monte-Carlo resampling (\code{"approximate"}).
    Alternatively, the functions 
    \code{\link{approximate}} or \code{\link{asymptotic}} can be
    used to specify how the exact conditional distribution of the test statistic
    should be calculated or approximated.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the independence of \code{y} and \code{x} is
  tested, \code{block} defines an optional factor for stratification. 
  \code{chisq_test} implements Pearson's chi-squared test, 
   \code{cmh_test} the Cochran-Mantel-Haenzsel
  test and \code{lbl_test} the linear-by-linear association test for ordered
  data.

  In case either \code{x} or \code{y} are ordered factors, the
  corresponding linear-by-linear association test is performed by all the
  procedures.
  \code{lbl_test} coerces factors to class \code{ordered} under any
  circumstances. The default scores are \code{1:nlevels(x)} and 
  \code{1:nlevels(y)}, respectively. The default scores can be changed 
  via the \code{scores} argument (see \code{\link{independence_test}}), 
  for example
  \code{scores = list(y = 1:3, x = c(1, 4, 6))} first triggers a coercion
  to class \code{ordered} of both variables and attaches the list elements
  as scores to the corresponding factors. The length of a score vector needs 
  to be equal the number of levels of the factor of interest.

  The authoritative source for details on the documented test procedures
  is Agresti (2002).

}
\value{

  An object inheriting from class \code{\link{IndependenceTest-class}} with
  methods \code{\link{show}}, \code{\link{statistic}}, \code{\link{expectation}}, 
  \code{\link{covariance}} and \code{\link{pvalue}}. The null distribution 
  can be inspected by \code{\link{pperm}}, \code{\link{dperm}},
  \code{\link{qperm}} and \code{\link{support}} methods.

}
\references{

    Alan Agresti (2002), \emph{Categorical Data Analysis}. Hoboken, New
    Jersey: John Wiley & Sons.

}
\examples{

\dontshow{
    set.seed(290875)
}

data("jobsatisfaction", package = "coin")

### for females only
chisq_test(as.table(jobsatisfaction[,,"Female"]), 
    distribution = approximate(B = 9999))

### both Income and Job.Satisfaction unordered
cmh_test(jobsatisfaction)

### both Income and Job.Satisfaction ordered, default scores
lbl_test(jobsatisfaction)

### both Income and Job.Satisfaction ordered, alternative scores
lbl_test(jobsatisfaction, scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                                        Income = c(3, 10, 20, 35)))

### the same, null distribution approximated
cmh_test(jobsatisfaction, scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                                        Income = c(3, 10, 20, 35)),
         distribution = approximate(B = 10000))

### Smoking and HDL cholesterin status
### (from Jeong, Jhun and Kim, 2005, CSDA 48, 623-631, Table 2)
smokingHDL <- as.table(
    matrix(c(15,  8, 11,  5, 
              3,  4,  6,  1, 
              6,  7, 15, 11, 
              1,  2,  3,  5), ncol = 4,
           dimnames = list(smoking = c("none", "< 5", "< 10", ">=10"), 
                           HDL = c("normal", "low", "borderline", "abnormal"))
    ))
### use interval mid-points as scores for smoking
lbl_test(smokingHDL, scores = list(smoking = c(0, 2.5, 7.5, 15)))

### Cochran-Armitage trend test for proportions
### Lung tumors in female mice exposed to 1,2-dichloroethane
### Encyclopedia of Biostatistics (Armitage & Colton, 1998), 
### Chapter Trend Test for Counts and Proportions, page 4578, Table 2
lungtumor <- data.frame(dose = rep(c(0, 1, 2), c(40, 50, 48)),
                        tumor = c(rep(c(0, 1), c(38, 2)),
                                  rep(c(0, 1), c(43, 7)),
                                  rep(c(0, 1), c(33, 15))))
table(lungtumor$dose, lungtumor$tumor)

### Cochran-Armitage test (permutation equivalent to correlation 
### between dose and tumor), cf. Table 2 for results
independence_test(tumor ~ dose, data = lungtumor, teststat = "quad")

### linear-by-linear association test with scores 0, 1, 2
### is identical with Cochran-Armitage test
lungtumor$dose <- ordered(lungtumor$dose)
independence_test(tumor ~ dose, data = lungtumor, teststat = "quad",
                  scores = list(dose = c(0, 1, 2)))

}
\keyword{htest}
