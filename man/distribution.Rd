\name{distribution}
\alias{exact}
\alias{asymptotic}
\alias{approximate}
\title{ Distribution under the Null Hypothesis }
\description{

    Specification of the exact, or approximation of the exact,
    conditional distribution of test statistics under the null hypothesis.

}
\usage{
exact(algorithm = c("shift", "split-up"), fact = NULL)
approximate(B = 1000)
asymptotic(maxpts = 25000, abseps = 0.001, releps = 0)
}
\arguments{
  \item{algorithm}{a character, specifying the algorithm to be used
                   for the computation of the exact conditional distribution.}
  \item{fact}{a positive integer to multiply the response values with. }
  \item{B}{a positive integer, the number of Monte-Carlo replications
           to approximate the exact conditional distribution.}
  \item{maxpts}{a positive integer, the maximum number of function 
                values, see \code{\link[mvtnorm]{pmvnorm}}.}
  \item{abseps}{double, the absolute error tolerance, 
                see \code{\link[mvtnorm]{pmvnorm}}.}
  \item{releps}{double, relative error tolerance, 
                see \code{\link[mvtnorm]{pmvnorm}}.}  

}
\details{

  The \code{distribution} argument to \code{\link{independence_test}} 
  can be specified with additional arguments using those functions.

}
\value{
  An object of class \code{exact}, \code{approximate} or \code{asymptotic},
  respectively.
}
\keyword{htest}
