\name{PermutationDistribution-methods}
\docType{methods}
\alias{dperm}
\alias{dperm-methods}
\alias{dperm,NullDistribution-method}
\alias{dperm,IndependenceTest-method}
\alias{pperm}
\alias{pperm-methods}
\alias{pperm,NullDistribution-method}
\alias{pperm,IndependenceTest-method}
\alias{qperm}
\alias{qperm-methods}
\alias{qperm,NullDistribution-method}
\alias{qperm,IndependenceTest-method}
\alias{rperm}
\alias{rperm-methods}
\alias{rperm,NullDistribution-method}
\alias{rperm,IndependenceTest-method}
\alias{support}
\alias{support-methods}
\alias{support,NullDistribution-method}
\alias{support,IndependenceTest-method}
\title{Computation of the Permutation Distribution}
\description{
  Methods for computation of the density function, distribution function,
  quantile function, random numbers and support of the permutation distribution.
}
\usage{
\S4method{dperm}{NullDistribution}(object, x, ...)
\S4method{dperm}{IndependenceTest}(object, x, ...)

\S4method{pperm}{NullDistribution}(object, q, ...)
\S4method{pperm}{IndependenceTest}(object, q, ...)

\S4method{qperm}{NullDistribution}(object, p, ...)
\S4method{qperm}{IndependenceTest}(object, p, ...)

\S4method{rperm}{NullDistribution}(object, n, ...)
\S4method{rperm}{IndependenceTest}(object, n, ...)

\S4method{support}{NullDistribution}(object, ...)
\S4method{support}{IndependenceTest}(object, ...)
}
\arguments{
  \item{object}{
    an object from which the density function, distribution function, quantile
    function, random numbers or support of the permutation distribution can be
    computed.
  }
  \item{x, q}{
    a numeric vector, the quantiles for which the density function or
    distribution function is computed.
  }
  \item{p}{
    a numeric vector, the probabilities for which the quantile function is
    computed.
  }
  \item{n}{
    a numeric vector, the number of observations.  If \code{length(n) > 1}, the
    length is taken to be the number required.
  }
  \item{...}{
    further arguments to be passed to methods.
  }
}
\details{
  The methods \code{dperm}, \code{pperm}, \code{qperm}, \code{rperm} and
  \code{support} compute the density function, distribution function, quantile
  function, random deviates and support, respectively, of the permutation
  distribution.
}
\value{
  The density function, distribution function, quantile function, random
  deviates or support of the permutation distribution computed from
  \code{object}.  A numeric vector.
}
\note{
  The density of asymptotic permutation distributions for maximum-type tests or
  exact permutation distributions obtained by the split-up algorithm is reported
  as \code{NA}.  The quantile function of asymptotic permutation distributions
  for maximum-type tests cannot be computed for \code{p} less than 0.5, due to
  limitations in the \pkg{mvtnorm} package.  The support of exact permutation
  distributions obtained by the split-up algorithm is reported as \code{NA}.

  In versions of \pkg{coin} prior to 1.1-0, the support of asymptotic
  permutation distributions was given as an interval containing 99.999 \% of the
  probability mass.  It is now reported as \code{NA}.
}
\examples{
## Two-sample problem
dta <- data.frame(
    y = rnorm(20),
    x = gl(2, 10)
)

## Exact Ansari-Bradley test
at <- ansari_test(y ~ x, data = dta, distribution = "exact")

## Support of the exact distribution of the Ansari-Bradley statistic
supp <- support(at)

## Density of the exact distribution of the Ansari-Bradley statistic
dens <- dperm(at, x = supp)

## Plotting the density
plot(supp, dens, type = "s")

## 95\% quantile
qperm(at, p = 0.95)

## One-sided p-value
pperm(at, q = statistic(at))

## Random number generation
rperm(at, n = 5)
}
\keyword{methods}
\keyword{htest}
\keyword{distribution}
