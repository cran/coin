\name{PermutationDistribution-methods}
\docType{methods}
\alias{dperm}
\alias{dperm-methods}
\alias{dperm,AsymptNullDistribution-method}
\alias{dperm,IndependenceTest-method}
\alias{dperm,NullDistribution-method}
\alias{pperm}
\alias{pperm-methods}
\alias{pperm,AsymptNullDistribution-method}
\alias{pperm,IndependenceTest-method}
\alias{pperm,NullDistribution-method}
\alias{qperm}
\alias{qperm-methods}
\alias{qperm,AsymptNullDistribution-method}
\alias{qperm,IndependenceTest-method}
\alias{qperm,NullDistribution-method}
\alias{rperm}
\alias{rperm-methods}
\alias{rperm,IndependenceTest-method}
\alias{rperm,NullDistribution-method}
\alias{support}
\alias{support-methods}
\alias{support,IndependenceTest-method}
\alias{support,NullDistribution-method}
\title{Computation of the Permutation Distribution}
\description{
  Methods for computation of the density function, distribution function,
  quantile function, random numbers and support of the permutation distribution.
}
\usage{
\S4method{dperm}{AsymptNullDistribution}(object, x, \dots)
\S4method{dperm}{IndependenceTest}(object, x, \dots)
\S4method{dperm}{NullDistribution}(object, x, \dots)

\S4method{pperm}{AsymptNullDistribution}(object, q, \dots)
\S4method{pperm}{IndependenceTest}(object, q, \dots)
\S4method{pperm}{NullDistribution}(object, q, \dots)

\S4method{qperm}{AsymptNullDistribution}(object, p, \dots)
\S4method{qperm}{IndependenceTest}(object, p, \dots)
\S4method{qperm}{NullDistribution}(object, p, \dots)

\S4method{rperm}{IndependenceTest}(object, n, \dots)
\S4method{rperm}{NullDistribution}(object, n, \dots)

\S4method{support}{IndependenceTest}(object, \dots)
\S4method{support}{NullDistribution}(object, \dots)
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
  \item{\dots}{
    further arguments to be passed to methods.
  }
}
\details{
  The methods \code{dperm}, \code{pperm}, \code{qperm}, \code{rperm} and
  \code{support} compute the density function, distribution function, quantile
  function, random numbers and support, respectively, of the permutation
  distribution.
}
\value{
  The density function, distribution function, quantile function, random numbers
  or support of the permutation distribution computed from \code{object}.  A
  numeric vector.
}
\note{
  The density of asymptotic permutation distributions for maximum-type tests or
  exact permutation distributions obtained by the split-up algoritm is reported
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
dens <- dperm(at, supp)

## Plotting the density
plot(supp, dens, type = "s")

## 95 \% quantile
qperm(at, 0.95)

## One-sided p-value
pperm(at, statistic(at))

## Random number generation
rperm(at, 5)
}
\keyword{methods}
\keyword{htest}
\keyword{distribution}
