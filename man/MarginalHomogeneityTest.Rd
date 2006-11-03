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
\method{mh_test}{SymmetryProblem}(object, distribution = c("asymptotic", "approximate"), ...) 
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
    can be approximated by its asymptotic distribution (\code{asymptotic}) 
    or via Monte-Carlo resampling (\code{approximate}).
    Alternatively, the functions 
    \code{\link{approximate}} or \code{\link{asymptotic}} can be
    used to specify how the exact conditional distribution of the test statistic
    should be calculated or approximated.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of independence of row and column totals is tested.
  The corresponding test for binary factors \code{x} and \code{y} is known
  as McNemar test. For larger tables, Stuart's \eqn{W_0} 
  statistic (Stuart, 1955, Agresti, 2002, page 422, also known as Stuart-Maxwell test) 
  is computed. The marginal homogeneity
  statistic \eqn{W} of Bhapkar (1966) can be derived from \eqn{W_0} 
  via \eqn{W = W_0 / (1 - W_0 / n)} (see Agresti, 2002, page 422).

  Scores must be a list of length one (row and column scores coincide). When
  scores are given or if \code{x} is ordered, the corresponding 
  linear association test is computed (see Agresti, 2002).

  Note that for a large number of observations, this function 
  is rather inefficient.

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\references{

    Alan Agresti (2002). \emph{Categorical Data Analysis}. Hoboken, New  
    Jersey: John Wiley & Sons.

    V. P. Bhapkar (1966). A note on the equivalence of two test criteria for hypotheses
    in categorical data. \emph{Journal of the American Statistical Association} \bold{61}, 
    228--235.

    Alan Stuart (1955). A test for homogeneity of the marginal distributions in a two-way
    classification. \emph{Biometrika} \bold{42}(3/4), 412--416.

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

### example taken from 
### http://ourworld.compuserve.com/homepages/jsuebersax/mcnemar.htm
rating <- c("low", "moderate", "high")
x <- as.table(matrix(c(20, 10,  5,
                       3, 30, 15,
                       0,  5, 40), 
                     ncol = 3, byrow = TRUE,
                     dimnames = list(Rater1 = rating, Rater2 = rating)))
### test statistic W_0 = 13.76
mh_test(x)

}
\keyword{htest}
