\name{jobsatisfaction}
\alias{jobsatisfaction}
\docType{data}
\title{ Income and Job Satisfaction }
\description{
  Income and job satisfaction by gender.
}
\usage{data(jobsatisfaction)}
\format{
 A tree-dimensional contigency table with dimensions
 \description{
     \item{Income}{a factor at levels \code{<5000, 5000-15000, 15000-25000}
                   and \code{>25000}.}
     \item{Job.Satisfaction}{a factor at levels \code{Very Dissatisfied, A
                             Little Dissatisfied, Moderately Satisfied} and 
                             \code{Very Satisfied}.}
     \item{Gender}{a factor at levels \code{Female} and \code{Male}.}
 }
}
\details{

  The data are given in Table 7.8, page 288, in Agresti (2002).

}
\source{

    Alan Agresti (2002), \emph{Categorical Data Analysis}. Hoboken, New
    Jersey: John Wiley & Sons.

}
\examples{

data(jobsatisfaction, package = "coin")

# Generalized Cochran-Mantel-Haenzel test
cmh_test(jobsatisfaction)

}
\keyword{datasets}
