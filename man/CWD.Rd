\name{CWD}
\docType{data}
\alias{CWD}
\title{Coarse Woody Debris}
\description{
  Carbon flux on six pieces of wood.
}
\usage{CWD}
\format{
  A data frame with 13 observations on 8 variables.
  \describe{
    \item{\code{sample2}}{
      carbon flux measurement for 2nd piece of wood.
    }
    \item{\code{sample3}}{
      carbon flux measurement for 3rd piece of wood.
    }
    \item{\code{sample4}}{
      carbon flux measurement for 4th piece of wood.
    }
    \item{\code{sample6}}{
      carbon flux measurement for 6th piece of wood.
    }
    \item{\code{sample7}}{
      carbon flux measurement for 7th piece of wood.
    }
    \item{\code{sample8}}{
      carbon flux measurement for 8th piece of wood.
    }
    \item{\code{trend}}{
      measurement day (in days from beginning).
    }
    \item{\code{time}}{
      date of measurement.
    }
  }
}
\details{
  Coarse woody debris (CWD, dead wood greater than 10 cm in diameter) is
  a large stock of carbon in tropical forests, yet the flux of carbon
  out of this pool, via respiration, is poorly resolved
  \bibcitep{coin::chambers_2001}.  The heterotrophic process involved in CWD
  respiration should respond to reductions in moisture availability,
  which occurs during dry season \bibcitep{coin::chambers_2001}.

  CWD respiration measurements were taken in a tropical forest in west French
  Guiana, which experiences extreme contrasts in wet and dry season
  \bibcitep{coin::bonal_2008}. An infrared gas analyzer and a clear chamber sealed to
  the wood surface were used to measure the flux of carbon out of the wood
  \bibcitep{coin::stahl_2011}.  Measurements were repeated 13 times, from July
  to November 2011, on six pieces of wood during the transition into and out of
  the dry season.  The aim is to assess if there were shifts in the CWD
  respiration of any of the pieces in response to the transition into (early
  August) and out of (late October) the dry season.

  \bibcitet{coin::Zeileis_Hothorn_2013} investigated the six-variate series of CO\eqn{_2}
  reflux, aiming to find out whether the reflux had changed over the sampling
  period in at least one of the six wood pieces.
}
\references{\bibshow{*}}
\source{
  The coarse woody debris respiration data were kindly provided by Lucy Rowland
  (School of GeoSciences, University of Edinburgh).
}
\examples{
## Zeileis and Hothorn (2013, pp. 942-944)
## Approximative (Monte Carlo) maximally selected statistics
CWD[1:6] <- 100 * CWD[1:6] # scaling (to avoid harmless warning)
mt <- maxstat_test(sample2 + sample3 + sample4 +
                   sample6 + sample7 + sample8 ~ trend, data = CWD,
                   distribution = approximate(nresample = 100000))

## Absolute maximum of standardized statistics (t = 3.08)
statistic(mt)

## 5\% critical value (t_0.05 = 2.86)
(c <- qperm(mt, 0.95))

## Only 'sample8' exceeds the 5\% critical value
sts <- statistic(mt, type = "standardized")
idx <- which(sts > c, arr.ind = TRUE)
sts[unique(idx[, 1]), unique(idx[, 2]), drop = FALSE]
}
\keyword{datasets}
