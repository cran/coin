\name{MarginalHomogeneityTest}
\alias{mh_test}
\alias{mh_test.table}
\alias{mh_test.formula}
\alias{mh_test.SymmetryProblem}
\title{ Marginal Homogeneity Test }
\description{
    Testing marginal homogeneity in a complete block design.
}
\usage{
\method{mh_test}{formula}(formula, data, subset = NULL, \dots)
\method{mh_test}{table}(object, ...)
\method{mh_test}{SymmetryProblem}(object, distribution = c("asympt", "approx"), ...) 
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a factor giving the data values and 
    \code{x} a factor with two or more levels giving the corresponding
    replications. \code{block} is an
    optional factor (which is generated automatically when omitted).}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{object}{an object inheriting from class \code{SymmetryProblem} or a 
           \code{table} with identical \code{dimnames} attributes.}
  \item{distribution}{a character, the null distribution of the test statistic
    can be approximated by its asymptotic distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of independence of row and column totals is tested.
  The corresponding test for binary factors \code{x} and \code{y} is known
  as McNemar test.

  Scores must be a list of length one (row and column scores coincide). When
  scores are given or if \code{x} is ordered, the corresponding 
  linear association test is computed (see Agresti, 2002).

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\references{

    Alan Agresti (2002), \emph{Categorical Data Analysis}. Hoboken, New  
    Jersey: John Wiley & Sons.

}
\examples{

### Opinions on Pre- and Extramarital Sex, Agresti (2002), page 421
opinions <- c("always wrong", "almost always wrong", 
              "wrong only sometimes", "not wrong at all")

PreExSex <- as.table(matrix(c(144, 33, 84, 126, 
                                2,  4, 14,  29, 
                                0,  2,  6,  25, 
                                0,  0,  1,  5), nrow = 4, 
                            dimnames = list(PremaritalSex = opinions,
                                            ExtramaritalSex = opinions)))

### treating response as nominal
mh_test(PreExSex)

### and as ordinal
mh_test(PreExSex, scores = list(response = 1:length(opinions)))

}
\keyword{htest}
