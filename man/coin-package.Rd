\name{coin-package}
\docType{package}
\alias{coin-package}
\alias{coin}
\title{General Information on the \pkg{coin} Package}
\description{
  The \pkg{coin} package provides an implementation of a general framework for
  conditional inference procedures commonly known as \emph{permutation tests}.
  The framework was developed by \bibcitet{coin::strasserweber1999} and is based on a
  multivariate linear statistic and its conditional expectation, covariance and
  limiting distribution.  These results are utilized to construct tests of
  independence between two sets of variables.

  The package does not only provide a flexible implementation of the abstract
  framework, but also provides a large set of convenience functions implementing
  well-known as well as lesser-known classical and non-classical test procedures
  within the framework.  Many of the tests presented in prominent text books,
  such as \bibcitet{coin::hollanderwolfe1999} or
  \bibcitet{coin::agresti2002}, are immediately
  available or can be implemented without much effort.  Examples include linear
  rank statistics for the two- and \eqn{K}-sample location and scale problem
  against ordered and unordered alternatives including post-hoc tests for
  arbitrary contrasts, tests of independence for contingency tables, two- and
  \eqn{K}-sample tests for censored data, tests of independence between two
  continuous variables as well as tests of marginal homogeneity and symmetry.
  Approximations of the exact null distribution via the limiting distribution or
  conditional Monte Carlo resampling are available for every test procedure,
  while the exact null distribution is currently available for univariate
  two-sample problems only.

  The salient parts of the Strasser-Weber framework are elucidated by
  \bibcitet{coin::hothorn:2006:amstat} and a thorough description of the software implementation
  is given by \bibcitet{coin::hothorn+hornik+vandewiel:2008}.
}
\author{
  This package is authored by \cr
  Torsten Hothorn <Torsten.Hothorn@R-project.org>, \cr
  Kurt Hornik <Kurt.Hornik@R-project.org>, \cr
  Mark A. van de Wiel <Mark.vdWiel@vumc.nl>, \cr
  Henric Winell <Henric.Winell@statistics.uu.se> and \cr
  Achim Zeileis <Achim.Zeileis@R-project.org>.
}
\references{\bibshow{*}}
\examples{
\dontrun{
## Generate doxygen documentation if you are interested in the internals:
## Download source package into a temporary directory
tmpdir <- tempdir()
tgz <- download.packages("coin", destdir = tmpdir, type = "source")[2]
## Extract contents
untar(tgz, exdir = tmpdir)
## Run doxygen (assuming it is installed)
wd <- setwd(file.path(tmpdir, "coin"))
system("doxygen inst/doxygen.cfg")
setwd(wd)
## Have fun!
browseURL(file.path(tmpdir, "coin", "inst",
                    "documentation", "html", "index.html"))}
}
\keyword{package}
