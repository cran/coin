\name{jobsatisfaction}
\docType{data}
\alias{jobsatisfaction}
\title{Income and Job Satisfaction}
\description{
  Income and job satisfaction by gender.
}
\usage{jobsatisfaction}
\format{
  A contingency table with 104 observations on 3 variables.
  \describe{
    \item{\code{Income}}{
      a factor with levels \code{"<5000"}, \code{"5000-15000"},
      \code{"15000-25000"} and \code{">25000"}.
    }
    \item{\code{Job.Satisfaction}}{
      a factor with levels \code{"Very Dissatisfied"},
      \code{"A Little Satisfied"}, \code{"Moderately Satisfied"} and
      \code{"Very Satisfied"}.
    }
    \item{\code{Gender}}{
      a factor with levels \code{"Female"} and \code{"Male"}.
    }
  }
}
\details{
  This data set was given in \bibcitet{|coin::agresti2002|Tab. 7.8 on p. 288}.  
  \bibcitet{coin::winell_2018} used the data to demonstrate a
  score-independent test for ordered categorical data.
}
\references{\bibshow{*}}
\examples{
## Approximative (Monte Carlo) linear-by-linear association test
lbl_test(jobsatisfaction, distribution = approximate(nresample = 10000))

\dontrun{
## Approximative (Monte Carlo) score-independent test
## Winell and Lindbaeck (2018)
(it <- independence_test(jobsatisfaction,
                         distribution = approximate(nresample = 10000),
                         xtrafo = function(data)
                             trafo(data, factor_trafo = function(x)
                                 zheng_trafo(as.ordered(x))),
                         ytrafo = function(data)
                             trafo(data, factor_trafo = function(y)
                                 zheng_trafo(as.ordered(y)))))

## Extract the "best" set of scores
ss <- statistic(it, type = "standardized")
idx <- which(abs(ss) == max(abs(ss)), arr.ind = TRUE)
ss[idx[1], idx[2], drop = FALSE]}
}
\keyword{datasets}
