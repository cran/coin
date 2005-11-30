\name{mercuryfish}
\alias{mercuryfish}
\docType{data}
\title{ Chromosomal Effects of Mercury Contaminated Fish Consumption }
\description{
    The mercury level in the blood, the proportion of cells with
    abnormalities and the proportion of cells with chromosome aberrations
    for a group of consuments of mercury contaminated fish and a control
    group.
}
\usage{data("mercuryfish")}
\format{
  A data frame with 39 observations on the following 4 variables.
  \describe{
    \item{group}{a factor with levels \code{control} and \code{exposed}.}
    \item{mercury}{the level of mercury in the blood.}
    \item{abnormal}{the proportion of cells with structural abnormalities.}
    \item{ccells}{the proportion of cells with asymmetrical or
                  incomplete-symmetrical chromosome aberrations called 
                  \eqn{C_u} cells.}
  }
}
\details{

  Subjects who ate contaminated fish for more than three years in the
  \code{exposed} group and subjects of a control group are to be compared.
  Instead of a multivariate comparison, Rosenbaum (1994) 
  applied a coherence criterion. The observations are partially ordered: an
  observation is smaller than another when all three variables (\code{mercury},
  \code{abnormal} and \code{ccells}) are smaller and a score reflecting the
  `ranking' is attached to each observation. The distribution of the scores
  in both groups is to be compared and the corresponding test is called
  `POSET-test' (partially ordered sets).
  
}
\source{

    S. Skerfving, K. Hansson, C. Mangs, J. Lindsten, N. Ryman (1974),
    Methylmercury-induced chromosome damage in men.
    \emph{Environmental Research} \bold{7}, 83--98.
}
\references{

    P. R. Rosenbaum (1994), Coherence in Observational Studies.
    \emph{Biometrics} \bold{50}, 368--374.

}
\examples{
data("mercuryfish")

coherence <- function(data) {
    x <- as.matrix(data)
    matrix(apply(x, 1, function(y)
        sum(colSums(t(x) < y) == ncol(x)) - 
        sum(colSums(t(x) > y) == ncol(x))), ncol = 1)
}

### POSET-test
poset <- independence_test(mercury + abnormal + ccells ~ group, data =
                           mercuryfish, ytrafo = coherence)

### linear statistic (T in Rosenbaum's, 1994, notation)
statistic(poset, "linear")

### expectation
expectation(poset)

### variance (there is a typo in Rosenbaum, 1994, page 371, 
### last paragraph Section 2)
covariance(poset)

### the standardized statistic
statistic(poset)

### and asymptotic p-value
pvalue(poset)

### exact p-value
independence_test(mercury + abnormal + ccells ~ group, data =
                  mercuryfish, ytrafo = coherence, distribution = "exact")

### multivariate analysis
mvtest <- independence_test(mercury + abnormal + ccells ~ group, 
                            data = mercuryfish)

### global p-value
pvalue(mvtest)

### adjusted univariate p-value
pvalue(mvtest, method = "single-step")

}
\keyword{datasets}
