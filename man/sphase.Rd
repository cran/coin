\name{sphase}
\alias{sphase}
\non_function{}
\title{ S-phase Fraction of Tumor Cells }
\usage{data("sphase")}
\description{
  S-phase fraction of tumor cells in breast cancer patients. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{SPF}{S-phase fraction}
    \item{RFS}{recurrence free survival}
    \item{event}{censoring indicator: \code{FALSE} means
                 censored, \code{TRUE} is an event.}
  }
}
\details{

  The data have been used to address the question whether a simple cutpoint 
  in S-phase fraction can be used to discriminate between patients with
  good and bad prognosis (for example in Hothorn & Lausen, 2003).

}
\source{


    J. Pfisterer, F. Kommoss, W. Sauerbrei, D. Menzel, M. Kiechle, E. Giese,
    M. Hilgarth \& A. Pfleiderer (1995). DNA flow cytometry in node positive
    breast cancer: Prognostic value and correlation to morphological and
    clinical factors. \emph{Analytical and Quantitative Cytology and Histology}
    \bold{7}(6), 406--412.

}
\references{

   Torsten Hothorn & Berthold Lausen (2003).
   On the Exact Distribution of Maximally Selected Rank
   Statistics. \emph{Computational Statistics \& Data Analysis}
   \bold{43}, 121--137.

}
\examples{

  maxstat_test(Surv(RFS, event) ~ SPF, data = sphase)

  ### reproduce the test statistic reported in Hothorn & Lausen (2003)
  maxstat_test(Surv(RFS, event) ~ SPF, data = sphase, 
      ytrafo = function(data) trafo(data, surv_trafo = function(x) 
          logrank_trafo(x, ties.method = "HL")))

}
\keyword{datasets}
