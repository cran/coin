\name{ContingencyTests}
\alias{chisq_test}
\alias{chisq_test.formula}
\alias{chisq_test.table}
\alias{chisq_test.IndependenceProblem}
\alias{cmh_test}
\alias{cmh_test.formula}
\alias{cmh_test.table}
\alias{cmh_test.IndependenceProblem}
\alias{lbl_test}
\alias{lbl_test.formula}
\alias{lbl_test.table}
\alias{lbl_test.IndependenceProblem}
\concept{Pearson chi-squared test}
\concept{Generalized Cochran-Mantel-Haenszel test}
\concept{Linear-by-linear association test}
\concept{Cochran-Armitage test}
\title{Tests of Independence in Two- or Three-Way Contingency Tables}
\description{
  Testing the independence of two nominal or ordered factors.
}
\usage{
\method{chisq_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{chisq_test}{table}(object, ...)
\method{chisq_test}{IndependenceProblem}(object, ...)

\method{cmh_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{cmh_test}{table}(object, ...)
\method{cmh_test}{IndependenceProblem}(object, ...)

\method{lbl_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{lbl_test}{table}(object, ...)
\method{lbl_test}{IndependenceProblem}(object, ...)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y ~ x | block} where \code{y} and \code{x} are
    factors and \code{block} is an optional factor for stratification.
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
    an object inheriting from classes \code{"table"} or
    \code{"\linkS4class{IndependenceProblem}"}.
  }
  \item{...}{
    further arguments to be passed to \code{\link{independence_test}()}.
  }
}
\details{
  \code{chisq_test()}, \code{cmh_test()} and \code{lbl_test()} provide the
  Pearson chi-squared test, the generalized Cochran-Mantel-Haenszel test and the
  linear-by-linear association test.  A general description of these methods is
  given by \bibcitet{coin::agresti2002}.

  The null hypothesis of independence, or conditional independence given
  \code{block}, between \code{y} and \code{x} is tested.

  If \code{y} and/or \code{x} are ordered factors, the default scores,
  \code{1:nlevels(y)} and \code{1:nlevels(x)}, respectively, can be altered
  using the \code{scores} argument (see \code{\link{independence_test}()}); this
  argument can also be used to coerce nominal factors to class \code{"ordered"}.
  (\code{lbl_test()} coerces to class \code{"ordered"} under any circumstances.)
  If both \code{y} and \code{x} are ordered factors, a linear-by-linear
  association test is computed and the direction of the alternative hypothesis
  can be specified using the \code{alternative} argument.  For the Pearson
  chi-squared test, this extension was given by \bibcitet{coin::yates_1948} who also discussed
  the situation when either the response or the covariate is an ordered factor;
  see also \bibcitet{coin::cochran_1954} and \bibcitet{coin::armitage:1955} for the particular case when
  \code{y} is a binary factor and \code{x} is ordered.  The Mantel-Haenszel
  statistic \bibcitep{coin::Mantel+Haenszel:1959} was similarly extended by
  \bibcitet{coin::mantel_1963} and \bibcitet{coin::landis_1978}.

  The conditional null distribution of the test statistic is used to obtain
  \eqn{p}-values and an asymptotic approximation of the exact distribution is
  used by default (\code{distribution = "asymptotic"}).  Alternatively, the
  distribution can be approximated via Monte Carlo resampling or computed
  exactly for univariate two-sample problems by setting \code{distribution} to
  \code{"approximate"} or \code{"exact"}, respectively.  See
  \code{\link{asymptotic}()}, \code{\link{approximate}()} and
  \code{\link{exact}()} for details.
}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\note{
  The exact versions of the Pearson chi-squared test and the generalized
  Cochran-Mantel-Haenszel test do not necessarily result in the same
  \eqn{p}-value as Fisher's exact test \bibcitep{coin::davis_1986}.
}
\references{\bibshow{*}}
\examples{
## Example data
## Davis (1986, p. 140)
davis <- matrix(
    c(3,  6,
      2, 19),
    nrow = 2, byrow = TRUE
)
davis <- as.table(davis)

## Asymptotic Pearson chi-squared test
chisq_test(davis)
chisq.test(davis, correct = FALSE) # same as above

## Approximative (Monte Carlo) Pearson chi-squared test
ct <- chisq_test(davis,
                 distribution = approximate(nresample = 10000))
pvalue(ct)             # standard p-value
midpvalue(ct)          # mid-p-value
pvalue_interval(ct)    # p-value interval
size(ct, alpha = 0.05) # test size at alpha = 0.05 using the p-value

## Exact Pearson chi-squared test (Davis, 1986)
## Note: disagrees with Fisher's exact test
ct <- chisq_test(davis,
                 distribution = "exact")
pvalue(ct)             # standard p-value
midpvalue(ct)          # mid-p-value
pvalue_interval(ct)    # p-value interval
size(ct, alpha = 0.05) # test size at alpha = 0.05 using the p-value
fisher.test(davis)


## Laryngeal cancer data
## Agresti (2002, p. 107, Tab. 3.13)
cancer <- matrix(
    c(21, 2,
      15, 3),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        "Treatment" = c("Surgery", "Radiation"),
           "Cancer" = c("Controlled", "Not Controlled")
    )
)
cancer <- as.table(cancer)

## Exact Pearson chi-squared test (Agresti, 2002, p. 108, Tab. 3.14)
## Note: agrees with Fishers's exact test
(ct <- chisq_test(cancer,
                  distribution = "exact"))
midpvalue(ct)          # mid-p-value
pvalue_interval(ct)    # p-value interval
size(ct, alpha = 0.05) # test size at alpha = 0.05 using the p-value
fisher.test(cancer)


## Homework conditions and teacher's rating
## Yates (1948, Tab. 1)
yates <- matrix(
    c(141, 67, 114, 79, 39,
      131, 66, 143, 72, 35,
       36, 14,  38, 28, 16),
    byrow = TRUE, ncol = 5,
    dimnames = list(
           "Rating" = c("A", "B", "C"),
        "Condition" = c("A", "B", "C", "D", "E")
    )
)
yates <- as.table(yates)

## Asymptotic Pearson chi-squared test (Yates, 1948, p. 176)
chisq_test(yates)

## Asymptotic Pearson-Yates chi-squared test (Yates, 1948, pp. 180-181)
## Note: 'Rating' and 'Condition' as ordinal
(ct <- chisq_test(yates,
                  alternative = "less",
                  scores = list("Rating" = c(-1, 0, 1),
                                "Condition" = c(2, 1, 0, -1, -2))))
statistic(ct)^2 # chi^2 = 2.332

## Asymptotic Pearson-Yates chi-squared test (Yates, 1948, p. 181)
## Note: 'Rating' as ordinal
chisq_test(yates,
           scores = list("Rating" = c(-1, 0, 1))) # Q = 3.825


## Change in clinical condition and degree of infiltration
## Cochran (1954, Tab. 6)
cochran <- matrix(
    c(11,  7,
      27, 15,
      42, 16,
      53, 13,
      11,  1),
    byrow = TRUE, ncol = 2,
    dimnames = list(
              "Change" = c("Marked", "Moderate", "Slight",
                           "Stationary", "Worse"),
        "Infiltration" = c("0-7", "8-15")
    )
)
cochran <- as.table(cochran)

## Asymptotic Pearson chi-squared test (Cochran, 1954, p. 435)
chisq_test(cochran) # X^2 = 6.88

## Asymptotic Cochran-Armitage test (Cochran, 1954, p. 436)
## Note: 'Change' as ordinal
(ct <- chisq_test(cochran,
                  scores = list("Change" = c(3, 2, 1, 0, -1))))
statistic(ct)^2 # X^2 = 6.66


## Change in size of ulcer crater for two treatment groups
## Armitage (1955, Tab. 2)
armitage <- matrix(
    c( 6, 4, 10, 12,
      11, 8,  8,  5),
    byrow = TRUE, ncol = 4,
    dimnames = list(
        "Treatment" = c("A", "B"),
           "Crater" = c("Larger", "< 2/3 healed",
                        ">= 2/3 healed", "Healed")
    )
)
armitage <- as.table(armitage)

## Approximative (Monte Carlo) Pearson chi-squared test (Armitage, 1955, p. 379)
chisq_test(armitage,
           distribution = approximate(nresample = 10000)) # chi^2 = 5.91

## Approximative (Monte Carlo) Cochran-Armitage test (Armitage, 1955, p. 379)
(ct <- chisq_test(armitage,
                  distribution = approximate(nresample = 10000),
                  scores = list("Crater" = c(-1.5, -0.5, 0.5, 1.5))))
statistic(ct)^2 # chi_0^2 = 5.26


## Relationship between job satisfaction and income stratified by gender
## Agresti (2002, p. 288, Tab. 7.8)

## Asymptotic generalized Cochran-Mantel-Haenszel test (Agresti, p. 297)
(ct <- cmh_test(jobsatisfaction)) # CMH = 10.2001

## The standardized linear statistic
statistic(ct, type = "standardized")

## The standardized linear statistic for each block
statistic(ct, type = "standardized", partial = TRUE)

## Asymptotic generalized Cochran-Mantel-Haenszel test (Agresti, p. 297)
## Note: 'Job.Satisfaction' as ordinal
cmh_test(jobsatisfaction,
         scores = list("Job.Satisfaction" = c(1, 3, 4, 5))) # L^2 = 9.0342

## Asymptotic linear-by-linear association test (Agresti, p. 297)
## Note: 'Job.Satisfaction' and 'Income' as ordinal
(lt <- lbl_test(jobsatisfaction,
                scores = list("Job.Satisfaction" = c(1, 3, 4, 5),
                              "Income" = c(3, 10, 20, 35))))
statistic(lt)^2 # M^2 = 6.1563

## The standardized linear statistic
statistic(lt, type = "standardized")

## The standardized linear statistic for each block
statistic(lt, type = "standardized", partial = TRUE)
}
\keyword{htest}
