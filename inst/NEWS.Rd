\name{NEWS}
\encoding{UTF-8}
\title{News for Package \pkg{coin}}
\section{Changes in version 1.0-24 (2014-09-17)}{
  \itemize{
    \item as.integer(round(x)) is safer than as.integer(x) = trunc(x)
  }
}
\section{Changes in version 1.0-23 (2013-09-03)}{
  \itemize{
    \item fix \file{DESCRIPTION} and \file{NAMESPACE} issues
  }
}
\section{Changes in version 1.0-22 (2013-04-16)}{
  \itemize{
    \item documentation updates

    \item test update (for Sparc Solaris)
  }
}
\section{Changes in version 1.0-21 (2012-03-08)}{
  \itemize{
    \item move \file{inst/doc/vignettes/*} to \file{vignettes/*}

    \item new \code{CWD} data
  }
}
\section{Changes in version 1.0-20 (2011-10-04)}{
  \itemize{
    \item report a warning for rank tests with weights
  }
}
\section{Changes in version 1.0-19 (2011-08-17)}{
  \itemize{
    \item no \code{require()} in \code{.onLoad()}
  }
}
\section{Changes in version 1.0-18 (2010-12-21)}{
  \itemize{
    \item \code{exact()} can not work for \code{symmetry_test()}; spotted by LE
    PAPE Gilles \email{lepape.gilles@neuf.fr}
  }
}
\section{Changes in version 1.0-17 (2010-09-13}{
  \itemize{
    \item add dependencies required for \file{tests/*}
  }
}
\section{Changes in version 1.0-16 (2010-08-23)}{
  \itemize{
    \item change \code{ties.method = c("HollanderWolfe", "Pratt")} to
    \code{zero.method = c("Pratt", "Wilcoxon")} in \code{wilcoxsign_test()}
    following a suggestion by Fritz Scholz.  A warning is fired to avoid
    misleading and nonreproducible results.

    \item use \code{quantile(..., type = 1)} for computing approximate quantiles
    (suggested by Fritz Scholz)
  }
}
\section{Changes in version 1.0-14 (2010-08-10}{
  \itemize{
    \item documentation \code{Biobase::expressionSet}

    \item further checks and better error messages for nonsense data in (paired)
    two-sample tests
  }
}
\section{Changes in version 1.0-12 (2010-07-03)}{
  \itemize{
    \item better error message for
    \code{wilcoxsign_test(c(1, 1, 1) ~ c(1, 1, 1))}
  }
}
\section{Changes in version 1.0-11 (2010-04-13)}{
  \itemize{
    \item bugfix in discrete MTP adjustment
  }
}
\section{Changes in version 1.0-10 (2010-01-27)}{
  \itemize{
    \item improve upon \code{qperm()} for van de Wiel algorithm again
  }
}
\section{Changes in version 1.0-9 (2009-12-30)}{
  \itemize{
    \item again precision problems in van de Wiel algorithm: \code{qperm()}
    suffered from too large tolerances of the inverted probability function.
  }
}
\section{Changes in version 1.0-8 (2009-10-30)}{
  \itemize{
    \item precision problems in van de Wiel algorithm fixed.  In rare cases,
    \eqn{P(T < t)} was returned as \eqn{p}-value instead of \eqn{P(T \le t)}.
    Two statistics are now considered equal (in all computations) when the
    difference is smaller than \code{sqrt(.Machine$double.eps)}
  }
}
\section{Changes in version 1.0-7 (2009-09-29)}{
  \itemize{
    \item \code{dperm()}, \code{pperm()} and \code{qperm()} are vectorized also
    for exact and approximate null distributions
  }
}
\section{Changes in version 1.0-6 (2009-09-08)}{
  \itemize{
    \item exact distribution for independent two-sample problems with only two
    observations was wrong, spotted by Fritz Scholz
    \email{fscholz@u.washington.edu}
  }
}
\section{Changes in version 1.0-5 (2009-06-04)}{
  \itemize{
    \item documentation updates
  }
}
\section{Changes in version 1.0-4 (2009-05-11)}{
  \itemize{
    \item evaluate all \code{formula}e in \code{xxx_test()}'s
    \code{parent.frame}
  }
}
\section{Changes in version 1.0-3 (2009-04-01)}{
  \itemize{
    \item add alternative ties handling to \code{wilcoxsign_test()} feature
    request by Fritz Scholz \email{fscholz@u.washington.edu}

    \item vignette update
  }
}
\section{Changes in version 1.0-2 (2009-03-16)}{
  \itemize{
    \item fix location confidence interval problem spotted by Fritz Scholz
    \email{fscholz@u.washington.edu}
  }
}
\section{Changes in version 1.0-1 (2009-01-09)}{
  \itemize{
    \item add average scores for logrank test

    \item fix Rd problems
  }
}
\section{Changes in version 1.0-0 (2008-11-06)}{
  \itemize{
    \item JSS paper \url{http://www.jstatsoft.org/v28/i08/} documents version
    1.0-0

    \item set default for \code{<IndependenceTest>@method} to
    \samp{General Independence Test}

    \item \code{covariance()} always returns a covariance matrix

    \item add \code{show()} method for \code{"IndependenceTest"} objects

    \item functions supplied via \code{distribution} does not need to have a
    class

    \item export all classes

    \item compute linear statistics, expectations and (co)variances when
    constructing <IndependenceLinearStatistic> instead of
    <IndependenceTestStatistic> and define methods for this class.

    \item \code{maxstat_test()} with integer weights gave \code{NA} or wrong
    cutpoints
  }
}
\section{Changes in version 0.6-9 (2008-05-05)}{
  \itemize{
    \item \code{any()} was misused in one place, spotted by Kasper Daniel Hansen
    \email{khansen@stat.Berkeley.EDU}

    \item code cosmetics by Johannes Huesing \email{johannes@huesing.name}
  }
}
\section{Changes in version 0.6-7 (2007-10-09)}{
  \itemize{
    \item remove non_function entries from Rd files

    \item disable \code{Biobase} example in vignette
  }
}
\section{Changes in version 0.6-6 (2007-08-13)}{
  \itemize{
    \item fix precision problems with exact \eqn{p}-values (spotted by Michael
    Fay \email{mfay@niaid.nih.gov}
  }
}
\section{Changes in version 0.6-5 (2007-07-23)}{
  \itemize{
    \item add new vignette on technical details

    \item update to new \pkg{mvtnorm} 0.8-0

    \item check for overflow errors in \file{StreitbergRoehmel.c} (thanks for
    Michael Fay \email{mfay@niaid.nih.gov} for spotting this)
  }
}
\section{Changes in version 0.6-4 (2007-06-28)}{
  \itemize{
    \item add new argument \code{ordered_trafo} to \code{trafo()} and deal with
    ordered factors in a more transparent way (via a new function
    \code{of_trafo()})
  }
}
\section{Changes in version 0.6-3 (2007-06-21)}{
  \itemize{
    \item new vignette \sQuote{MAXtest}
  }
}
\section{Changes in version 0.6-2 (2007-05-13)}{
  \itemize{
    \item print name of \code{x} variable and levels when \code{x} is a factor

    \item new class \code{"IndependenceLinearStatistic"}
  }
}
\section{Changes in version 0.6-1 (2007-04-25)}{
  \itemize{
    \item \code{LazyLoad: yes}
  }
}
\section{Changes in version 0.6-0 (2007-03-01)}{
  \itemize{
    \item rename some C source files & update doxygen documentation

    \item simplify some methods and the class structure

    \item clean up code and simplify wrapper functions

    \item documentation updates

    \item deal with factors containing only one level

    \item do not compute design matrices in \code{ModelEnvFormula}

    \item expand weights if \code{distribution = approximate()}

    \item assign names to user-specified transformation in case they are missing

    \item the \code{distribution} argument may now be a function with one
    argument allowing for user-specified distributions.

    \item \code{surv_test()} ignored the \code{alternative} argument
  }
}
\section{Changes in version 0.5-2 (2006-11-03)}{
  \itemize{
    \item \code{maxstat_trafo} did not always choose the correct maximal
    cutpoint

    \item formula evaluation without \code{data} argument was partially broken
    (thanks to Achim Zeileis for spotting this)
  }
}
\section{Changes in version 0.5-1 (2006-10-27)}{
  \itemize{
    \item improve formula interface for \code{wilcoxsign_test()}
  }
}
\section{Changes in version 0.5-0 (2006-10-17)}{
  \itemize{
    \item implement exact distribution for symmetry problems (especially
    \code{wilcoxsign_test()})

    \item add more checks on StatXact examples
  }
}
\section{Changes in version 0.4-14 (2006-09-18)}{
  \itemize{
    \item fix problems reported by valgrind
  }
}
\section{Changes in version 0.4-13 (2006-09-13)}{
  \itemize{
    \item \file{coin-Ex.R} generated by \R 2.4.0

    \item Enhances: Biobase
  }
}
\section{Changes in version 0.4-12 (2006-09-07)}{
  \itemize{
    \item add $(FLIBS) to \file{Makevars}

    \item include doxygen documentation for C functions

    \item documentation updates
  }
}
\section{Changes in version 0.4-10 (2006-08-25)}{
  \itemize{
    \item \code{maxstat_test(y ~ x)} is now able to deal with unordered \code{x}
    variables

    \item \code{maxstat_test()} returns estimates of both the selected variable
    and the cutpoint in multivariate situations (as a list)

    \item more checks for \code{maxstat_test()}

    \item add \pkg{codetools} checks to tests
  }
}
\section{Changes in version 0.4-9 (2006-07-24)}{
  \itemize{
    \item fix typo: Homegeneity

    \item add ASA copyright to \sQuote{LegoCondInf} vignette
  }
}
\section{Changes in version 0.4-8 (2006-06-23)}{
  \itemize{
    \item clarify that the Stuart-Maxwell test is computed by \code{mh_test()}
    and add a further example (thanks to Henric Nilsson
    \email{henric.nilsson@phadia.com} for pointing this out!)
  }
}
\section{Changes in version 0.4-7 (2006-05-16)}{
  \itemize{
    \item improve \code{maxstat_trafo()}: some potential cutpoints could have
    been overlooked in case of ties equal to the maxprob and minprob sample
    quantiles

    \item reproduce two examples from Hothorn and Lausen (2003) in
    \code{?maxstat_test} (and added \code{hohnloser} data set)

    \item add more regression tests on \code{maxstat_test()}
  }
}
\section{Changes in version 0.4-6 (2006-04-24)}{
  \itemize{
    \item updated \sQuote{LegoCondInf} vignette

    \item printed names of test statistics are now \samp{chi-squared},
    \samp{maxT}, or \samp{Z} instead of \samp{T}
  }
}
\section{Changes in version 0.4-5 (2006-04-12)}{
  \itemize{
    \item user-supplied transformations can return a vector which is coerced to
    \code{matrix(..., ncol = 1)} internally

    \item \code{ytrafo} and \code{xtrafo} may take functions like \code{rank()}
    directly.  However, the use of \code{trafo()} is recommended.

    \item argument \code{teststat} may take values \code{"max"}, \code{"quad"}
    or \code{"scalar"}

    \item internal reimplementation of score handling for ordered factors

    \item handling of multiple ordered factors / ordered factors in multivariate
    situations implemented

    \item \code{f_trafo()} returns a design matrix with \code{NA} rows in case
    missing values were present
  }
}
\section{Changes in version 0.4-4 (2006-02-03)}{
  \itemize{
    \item make R CMD check happy and move \file{src/README} to
    \file{inst/README}
  }
}
\section{Changes in version 0.4-3 (2006-01-12)}{
  \itemize{
    \item \code{isequal()} must not check equality of attributes
  }
}
\section{Changes in version 0.4-2 (2005-12-21)}{
  \itemize{
    \item new data set \code{alpha}

    \item data set \code{alzheimer} is now a \code{"data.frame"}, not a
    \code{"table"}

    \item new vignette with more applications
  }
}
\section{Changes in version 0.4-1 (2005-12-02)}{
  \itemize{
    \item documentation updates

    \item remove unused \code{setAs} definitions

    \item \code{expectation()}, \code{variance()} and \code{covariance()} return
    named vectors or matrices.  Names for objects returned by \code{statistic()}
    have been partially improved.

    \item all observations with \code{NA}s are removed now

    \item some internal improvements (removed unused code, avoid duplicated code
    chunks etc.)

    \item more checks on \code{NA}s, blocks and multiple ordinal variables (one
    univariate problems are currently allowed to have ordinal variables)
  }
}
\section{Changes in version 0.3-8 (2005-11-29)}{
  \itemize{
    \item export \code{independence_test.table()}

    \item add \code{conf.level} attribute for \code{MCp} objects when
    quasi-randomized Monte Carlo procedures (from \pkg{mvtnorm}) have been used

    \item add new \code{var_trafo} argument to \code{trafo()}

    \item documentation updates
  }
}
\section{Changes in version 0.3-7 (2005-11-23)}{
  \itemize{
    \item \code{statistic()} now returns correct linear and standardized
    statistics when scores are in play

    \item \code{alzheimer} data added

    \item \code{photocar} data added

    \item documentation updates

    \item \code{support()} and \code{dperm()} methods have been improved for
    asymptotic and approximate null distributions
  }
}
\section{Changes in version 0.3-6 (2005-11-17)}{
  \itemize{
    \item internal C function \code{R_MonteCarloIndependenceTest} returns a
    \eqn{pq \times B)} matrix instead of a list with B elements (the linear
    statistic for each random permutation of the data)

    \item various simplifications for the computations of adjusted
    \eqn{p}-values

    \item \code{I()} in \code{formula}e could still cause trouble with class
    \code{AsIs}, now fixed

    \item logical variables are now allowed in \code{formula}e (and are treated
    like \code{factor}s)

    \item fix \code{print} method for marginal homogeneity tests
    (\samp{...stratified by block...})

    \item \code{R_kronecker} in C now available

    \item C versions of \code{nrow()} and \code{ncol()} return 1 or length for
    vectors

    \item some internal optimisation

    \item C functions \code{nrow(x)} and \code{ncol(x)} return \code{LENGTH(x)}
    or \code{1} when \code{x} has no \code{dim} attribute

    \item new R-interface function \code{R_kronecker} to \code{C_kronecker}
    (which returns a _vector_!)

    \item speedup of variance computations in internal functions
  }
}
\section{Changes in version 0.3-3 (2005-09-07)}{
  \itemize{
    \item fix more problems reported by new \pkg{codetools}.  Try to work around
    the \code{terms(y ~ ., data = data.frame(1:10))} problem in \R 2.2.0.
  }
}
\section{Changes in version 0.3-2 (2005-08-31)}{
  \itemize{
    \item fix some problems reported by \pkg{codetools}
  }
}
\section{Changes in version 0.3-1 (2005-08-20)}{
  \itemize{
    \item \code{alternative = "less"} and \code{alternative = "greater"} are now
    defined for \code{"maxtype"} statistics as well

    \item one- and two-sided single-step and step-down max-\eqn{T} \eqn{p}-value
    adjustments are now available from the appropriate \code{pvalue()} method
    (NOTE: those procedures have not been tested carefully, yet)

    In addition, the Bonferroni-adjustment by Westfall & Wolfinger (1997) is
    available now.  Note that the interface to \code{pvalue()} changed slightly,
    \code{adjustment = TRUE} was replaced by \code{method = "single-step"}.

    \item more examples added to the vignette.
  }
}
\section{Changes in version 0.2-14 (2005-07-14)}{
  \itemize{
    \item \code{distribution = "approximate"} for \code{"maxtype"} statistics
    was wrong in case both \code{xtrafo} \emph{and} \code{ytrafo} were
    multivariate (the conditional expectation was computed incorrectly)
  }
}
\section{Changes in version 0.2-13 (2005-07-11)}{
  \itemize{
    \item \code{I()} in \code{formula}e returns objects of class \code{AsIs}
    which caused troubles in \code{trafo()}
  }
}
\section{Changes in version 0.2-12 (2005-06-15)}{
  \itemize{
    \item functionality for \code{formula} parsing and evaluation is now
    imported from package \pkg{modeltools}

    \item \code{show()} returns objects (of class \code{"htest"}, for example)
    invisibly, really

    \item \code{maxstat_trafo()} is much faster now and returns a matrix with
    both row- and column names set appropriately
  }
}
\section{Changes in version 0.2-11 (2005-06-02)}{
  \itemize{
    \item the \code{distribution} argument now takes the return values of
    functions \code{exact()}, \code{approximate()} or \code{asymptotic()} as
    well.  Those functions can be used to specify parameters, such as the number
    of Monte Carlo replications via
    \code{..., distribution = approximate(B = 9999), ...}

    \item \code{show()} returns objects (of class \code{"htest"}, for example)
    invisibly

    \item \code{expectation()} returns a vector, not a matrix

    \item new generic \code{variance()} for extracting the variance(s) of linear
    statistics

    \item only variances (instead of the whole covariance matrix) is computed
    when the distribution of maximum type test statistics is to be
    \emph{approximated}

    \item \code{data} may be an object of class \code{"exprSet"} (->
    \pkg{Biobase}), the vignette has an example

    \item \code{logrank_trafo()} (and \code{surv_test()}) now have a
    \code{ties.method} argument, see \code{?surv_test} for more information
  }
}
\section{Changes in version 0.2-10 (2005-05-04)}{
  \itemize{
    \item \samp{asymptotical} \eqn{\rightarrow} \samp{asymptotic} in
    \code{print} methods

    \item \code{mercuryfish} example added

    \item \code{<x,y>trafo} now can return matrices with number of columns
    different from the lhs and rhs of \code{formula}
  }
}
\section{Changes in version 0.2-9 (2005-03-15)}{
  \itemize{
    \item \code{mergesort} is already defined in
    \file{/usr/include/stdlib.h}:270 on some platforms

    \item \code{delay} is deprecated in \R 2.1.0
  }
}
