\name{treepipit}
\docType{data}
\alias{treepipit}
\encoding{UTF-8}
\title{Tree Pipits in Franconian Oak Forests}
\description{
  Data on the population density of tree pipits, \emph{Anthus trivialis}, in
  Franconian oak forests including variables describing the forest ecosystem.
}
\usage{treepipit}
\format{
  A data frame with 86 observations on 10 variables.
  \describe{
    \item{\code{counts}}{
      the number of tree pipits observed.
    }
    \item{\code{age}}{
      age of the overstorey oaks taken from forest data.
    }
    \item{\code{coverstorey}}{
      cover of canopy overstorey (\%).  The crown cover is described relative to
      a fully stocked stand.  Very dense overstorey with multiple crown cover
      could reach values greater than 100\%.
    }
    \item{\code{coverregen}}{
      cover of regeneration and shrubs (\%).
    }
    \item{\code{meanregen}}{
      mean height of regeneration and shrubs.
    }
    \item{\code{coniferous}}{
      coniferous trees (\% per hectare).
    }
    \item{\code{deadtree}}{
      number of dead trees (per hectare).
    }
    \item{\code{cbpiles}}{
      number of crowns and branch piles (per hectare).  All laying crowns and
      branch piles were counted.  These were induced by logging and the creation
      of wind breaks.
    }
    \item{\code{ivytree}}{
      number of ivied trees (per hectare).
    }
    \item{\code{fdist}}{
      distance to the forest edge.  The closest distance to the forest edge was
      measured from the centre of each grid.
    }
  }
}
\details{
  This study is based on fieldwork conducted in three lowland oak forests in the
  Franconian region of northern Bavaria close to Uffenheim, Germany.  Diurnal
  breeding birds were sampled five times, from March to June 2002, using a
  quantitative grid mapping.  Each grid was a one-hectare square.  In total, 86
  sample sites were established in 9 stands.  All individuals were counted in
  time intervals of 7 min/grid during slow walks along the middle of the grid
  with a stop in the centre.  Environmental factors were measured for each
  grid.
}
\references{
  \enc{Müller}{Mueller}, J. and Hothorn, T.  (2004).  Maximally selected
  two-sample statistics as a new tool for the identification and assessment of
  habitat factors with an application to breeding bird communities in oak
  forests.  \emph{European Journal of Forest Research} \bold{123}(3), 219--228.
  \doi{10.1007/s10342-004-0035-5}
}
\examples{
## Asymptotic maximally selected statistics
maxstat_test(counts ~ age + coverstorey + coverregen + meanregen +
                      coniferous + deadtree + cbpiles + ivytree,
             data = treepipit)
}
\keyword{datasets}
