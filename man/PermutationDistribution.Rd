\name{PermutationDistribution}
\docType{methods}
\alias{pperm}
\alias{pperm-methods}
\alias{pperm,NullDistribution-method}
\alias{pperm,IndependenceTest-method}
\alias{pperm,ScalarIndependenceTest-method}
\alias{pperm,MaxTypeIndependenceTest-method}
\alias{pperm,QuadTypeIndependenceTest-method}
\alias{qperm}
\alias{qperm-methods}
\alias{qperm,NullDistribution-method}
\alias{qperm,IndependenceTest-method}
\alias{qperm,ScalarIndependenceTest-method}
\alias{qperm,MaxTypeIndependenceTest-method}
\alias{qperm,QuadTypeIndependenceTest-method}
\alias{dperm}
\alias{dperm-methods}
\alias{dperm,NullDistribution-method}
\alias{dperm,IndependenceTest-method}
\alias{dperm,ScalarIndependenceTest-method}
\alias{dperm,MaxTypeIndependenceTest-method}
\alias{dperm,QuadTypeIndependenceTest-method}
\alias{support}
\alias{support-methods}
\alias{support,NullDistribution-method}
\alias{support,IndependenceTest-method}
\alias{support,ScalarIndependenceTest-method}
\alias{support,MaxTypeIndependenceTest-method}
\alias{support,QuadTypeIndependenceTest-method}

\title{ Permutation Distribution of Conditional Independence Tests }

\description{
  Density, distribution, quantile and support of permutation distributions
  of conditional independence tests.
}
\usage{
dperm(object, x, ...)
pperm(object, q, ...)
qperm(object, p, ...)
support(object, p, ...)
}
\arguments{
  \item{object}{an object inheriting from class
    \code{\link{IndependenceTest-class}}.}
  \item{x}{value of a standardized statistic.}
  \item{q}{value of a standardized statistic.}
  \item{p}{a probability.}
  \item{\dots}{further arguments (currently ignored).}
}
\section{Methods}{
  \describe{
    \item{dperm}{density function of a permutation distribution.}
    \item{pperm}{distribution function of a permutation distribution.}
    \item{qperm}{quantile function for a permutation distribution.}
    \item{support}{the support of a permutation distribution.}
  }
}
\examples{

### artificial 2-sample problem
df <- data.frame(y = rnorm(20), x = gl(2, 10))

### Ansari-Bradley test
at <- ansari_test(y ~ x, data = df, distribution = "exact")

### density of the exact distribution of the Ansari-Bradley statistic
dens <- sapply(support(at), dperm, object = at)

### 95\% quantile
qperm(at, 0.95)

### one-sided p-value
pperm(at, statistic(at))

}
\keyword{methods}
\keyword{htest}
