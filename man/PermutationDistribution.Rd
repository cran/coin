\name{PermutationDistribution}
\docType{methods}
\alias{pperm}
\alias{pperm-methods}
\alias{pperm,NullDistribution-method}
\alias{pperm,AsymptNullDistribution-method}
\alias{pperm,IndependenceTest-method}
\alias{qperm}
\alias{qperm-methods}
\alias{qperm,NullDistribution-method}
\alias{qperm,AsymptNullDistribution-method}
\alias{qperm,IndependenceTest-method}
\alias{dperm}
\alias{dperm-methods}
\alias{dperm,NullDistribution-method}
\alias{dperm,AsymptNullDistribution-method}
\alias{dperm,IndependenceTest-method}
\alias{support}
\alias{support-methods}
\alias{support,NullDistribution-method}
\alias{support,IndependenceTest-method}

\title{ Permutation Distribution of Conditional Independence Tests }

\description{
  Density, distribution, quantile and support of permutation distributions
  of conditional independence tests.
}
\usage{
dperm(object, x, ...)
pperm(object, q, ...)
qperm(object, p, ...)
support(object, ...)
}
\arguments{
  \item{object}{an object inheriting from class
    \code{\link{IndependenceTest-class}}.}
  \item{x}{vector of a standardized statistics.}
  \item{q}{vector of a standardized statistics.}
  \item{p}{vector of probabilities.}
  \item{\dots}{further arguments to be passed to methods.}
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

  ### plot density
  plot(support(at), dens, type = "s")

  ### 95\% quantile
  qperm(at, 0.95)

  ### one-sided p-value
  pperm(at, statistic(at))

}
\keyword{methods}
\keyword{htest}
