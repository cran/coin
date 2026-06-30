\name{MarginalHomogeneityTests}
\alias{mh_test}
\alias{mh_test.formula}
\alias{mh_test.table}
\alias{mh_test.SymmetryProblem}
\concept{McNemar test}
\concept{Cochran Q test}
\concept{Stuart(-Maxwell) test}
\concept{Madansky test of interchangeability}
\title{Marginal Homogeneity Tests}
\description{
  Testing the marginal homogeneity of a repeated measurements factor in a
  complete block design.
}
\usage{
\method{mh_test}{formula}(formula, data, subset = NULL, weights = NULL, ...)
\method{mh_test}{table}(object, ...)
\method{mh_test}{SymmetryProblem}(object, ...)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y ~ x | block} where \code{y} and \code{x} are
    factors and \code{block} is an optional factor (which is generated
    automatically if omitted).
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
    weight for all observations.  (Not yet implemented!)
  }
  \item{object}{
    an object inheriting from classes \code{"table"} (with identical
    \code{dimnames} components) or \code{"\linkS4class{SymmetryProblem}"}.
  }
  \item{...}{
    further arguments to be passed to \code{\link{symmetry_test}()}.
  }
}
\details{
  \code{mh_test()} provides the McNemar test, the Cochran \eqn{Q} test, the
  Stuart(-Maxwell) test and the Madansky test of interchangeability.  A general
  description of these methods is given by \bibcitet{coin::agresti2002}.

  The null hypothesis of marginal homogeneity is tested.  The response variable
  and the measurement conditions are given by \code{y} and \code{x},
  respectively, and \code{block} is a factor where each level corresponds to
  exactly one subject with repeated measurements.

  This procedure is known as the McNemar test \bibcitep{coin::mcnemar_1947} when both \code{y}
  and \code{x} are binary factors, as the Cochran \eqn{Q} test
  \bibcitep{coin::cochran_1954}
  when \code{y} is a binary factor and \code{x} is a factor with an arbitrary
  number of levels, as the Stuart(-Maxwell) test
  \bibcitep{coin::stuart_1955,maxwell_1970}
  when \code{y} is a factor with an arbitrary number of levels and \code{x} is a
  binary factor, and as the Madansky test of interchangeability
  \bibcitep{coin::madansky_1963},
  which implies marginal homogeneity, when both \code{y} and \code{x} are
  factors with an arbitrary number of levels.
  %% See also Bhapkar (1966) for a more powerful Wald-type test.  (As noted by
  %% Ireland, Ku and Kullback (1969), Bhapkar's statistic is given by
  %% \eqn{Q / (1 - Q / n)}, where \eqn{Q} is Stuart's statistic and \eqn{n} is
  %% the sample size).

  If \code{y} and/or \code{x} are ordered factors, the default scores,
  \code{1:nlevels(y)} and \code{1:nlevels(x)}, respectively, can be altered
  using the \code{scores} argument (see \code{\link{symmetry_test}()}); this
  argument can also be used to coerce nominal factors to class \code{"ordered"}.
  If both \code{y} and \code{x} are ordered factors, a linear-by-linear
  association test is computed and the direction of the alternative hypothesis
  can be specified using the \code{alternative} argument.  This extension was
  given by \bibcitet{coin::birch_1965} who also discussed the situation when either the
  response or the measurement condition is an ordered factor; see also
  \bibcitet{coin::white_1982}.

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
  This function is currently computationally inefficient for data with a large
  number of pairs or sets.
}
\references{\bibshow{*}}
\examples{
## Performance of prime minister
## Agresti (2002, p. 409)
performance <- matrix(
    c(794, 150,
       86, 570),
    nrow = 2, byrow = TRUE,
    dimnames = list(
         "First" = c("Approve", "Disprove"),
        "Second" = c("Approve", "Disprove")
    )
)
performance <- as.table(performance)
diag(performance) <- 0 # speed-up: only off-diagonal elements contribute

## Asymptotic McNemar Test
mh_test(performance)

## Exact McNemar Test
mh_test(performance, distribution = "exact")


## Effectiveness of different media for the growth of diphtheria
## Cochran (1950, Tab. 2)
cases <- c(4, 2, 3, 1, 59)
n <- sum(cases)
cochran <- data.frame(
    diphtheria = factor(
        unlist(rep(list(c(1, 1, 1, 1),
                        c(1, 1, 0, 1),
                        c(0, 1, 1, 1),
                        c(0, 1, 0, 1),
                        c(0, 0, 0, 0)),
                   cases))
    ),
    media = factor(rep(LETTERS[1:4], n)),
    case =  factor(rep(seq_len(n), each = 4))
)

## Asymptotic Cochran Q test (Cochran, 1950, p. 260)
mh_test(diphtheria ~ media | case, data = cochran) # Q = 8.05

## Approximative Cochran Q test
mt <- mh_test(diphtheria ~ media | case, data = cochran,
              distribution = approximate(nresample = 10000))
pvalue(mt)             # standard p-value
midpvalue(mt)          # mid-p-value
pvalue_interval(mt)    # p-value interval
size(mt, alpha = 0.05) # test size at alpha = 0.05 using the p-value


## Opinions on Pre- and Extramarital Sex
## Agresti (2002, p. 421)
opinions <- c("Always wrong", "Almost always wrong",
              "Wrong only sometimes", "Not wrong at all")
PreExSex <- matrix(
    c(144, 33, 84, 126,
        2,  4, 14,  29,
        0,  2,  6,  25,
        0,  0,  1,   5),
    nrow = 4,
    dimnames = list(
          "Premarital Sex" = opinions,
        "Extramarital Sex" = opinions
    )
)
PreExSex <- as.table(PreExSex)

## Asymptotic Stuart test
mh_test(PreExSex)

## Asymptotic Stuart-Birch test
## Note: response as ordinal
mh_test(PreExSex, scores = list(response = 1:length(opinions)))


## Vote intention
## Madansky (1963, pp. 107-108)
vote <- array(
    c(120, 1,  8, 2,   2,  1, 2, 1,  7,
        6, 2,  1, 1, 103,  5, 1, 4,  8,
       20, 3, 31, 1,   6, 30, 2, 1, 81),
    dim = c(3, 3, 3),
    dimnames = list(
          "July" = c("Republican", "Democratic", "Uncertain"),
        "August" = c("Republican", "Democratic", "Uncertain"),
          "June" = c("Republican", "Democratic", "Uncertain")
    )
)
vote <- as.table(vote)

## Asymptotic Madansky test (Q = 70.77)
mh_test(vote)


## Cross-over study
## http://www.nesug.org/proceedings/nesug00/st/st9005.pdf
dysmenorrhea <- array(
    c(6, 2, 1,  3, 1, 0,  1, 2, 1,
      4, 3, 0, 13, 3, 0,  8, 1, 1,
      5, 2, 2, 10, 1, 0, 14, 2, 0),
    dim = c(3, 3, 3),
    dimnames =  list(
          "Placebo" = c("None", "Moderate", "Complete"),
         "Low dose" = c("None", "Moderate", "Complete"),
        "High dose" = c("None", "Moderate", "Complete")
    )
)
dysmenorrhea <- as.table(dysmenorrhea)

## Asymptotic Madansky-Birch test (Q = 53.76)
## Note: response as ordinal
mh_test(dysmenorrhea, scores = list(response = 1:3))

## Asymptotic Madansky-Birch test (Q = 47.29)
## Note: response and measurement conditions as ordinal
mh_test(dysmenorrhea, scores = list(response = 1:3,
                                    conditions = 1:3))
}
\keyword{htest}
