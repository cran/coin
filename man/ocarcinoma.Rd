\name{ocarcinoma}
\alias{ocarcinoma}
\docType{data}
\title{ Ovarian Carcinoma }
\description{
  Survival times of 35 women suffering ovarian carcinoma at stadium II and
IIA.
}
\usage{data(ocarcinoma)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{time}{time in days.}
    \item{event}{censoring indicator: \code{FALSE} means 
                 censored, \code{TRUE} is an event.}
    \item{stadium}{a factor with levels \code{II} and \code{IIA}. }
  }
}
\details{
  Data from Fleming et al. (1980, 1984), reanalysed in Schumacher and
Schulgen (2002).
}
\source{
  
    Thomas R. Fleming, Judith R. O'Fallon, Peter C. O'Brien & David 
    P. Harrington (1980),
    Modified Kolmogorov-Smirnov test procedures with applications to 
    arbitrarily censored data. \emph{Biometrics}, \bold{36}, 607--625.

    Thomas R. Fleming, Stephanie J. Green & David P. Harrington (1984),
    Considerations of monitoring and evaluating treatment effects in clinical
    trials. \emph{Controlled Clinical Trials}, \bold{5}, 55--66.

}
\references{

    Martin Schumacher & Gabi Schulgen (2002), \emph{Methodik klinischer 
    Studien: methodische Grundlagen der Planung, Durchf\"uhrung und
    Auswertung.} Springer, Heidelberg.

}
\examples{

data(ocarcinoma, package = "coin")

### logrank test with exact two-sided p-value
lrt <- surv_test(Surv(time, event) ~ stadium, data = ocarcinoma,
                    distribution = "exact")

### the test statistic
statistic(lrt)

### p-value
pvalue(lrt)

}
\keyword{datasets}
