\name{NullDistribution-class}
\docType{class}
\alias{NullDistribution-class}

\title{Class "NullDistribution" }
\description{ Representing reference distributions }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NullDistribution", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{q}:}{quantile function }
    \item{\code{d}:}{density function }
    \item{\code{support}:}{support of the distribution }
    \item{\code{parameters}:}{list of additional parameters }
    \item{\code{pvalue}:}{p value function }
    \item{\code{p}:}{distribution function }
    \item{\code{name}:}{character, name of the distribution }
  }
}
\section{Extends}{
Class \code{"\linkS4class{PValue}"}, directly.
}
\keyword{classes}
