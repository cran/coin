
### generic method for extracting pvalues from objects
setGeneric("pvalue", function(object, ...)
    standardGeneric("pvalue"))

setMethod(f = "pvalue",
          signature = "NullDistribution",
          definition = function(object, q, ...) {
              object@pvalue(q)
          }
)

setMethod(f = "pvalue",
          signature = "IndependenceTest",
          definition = function(object, q, ...) {
              pvalue(object@distribution, q = q)
          }
)

setMethod(f = "pvalue",
          signature = "ScalarIndependenceTest",
          definition = function(object, ...) {
              object@distribution@pvalue(object@statistic@teststatistic)
          }
)

setMethod(f = "pvalue",
          signature = "MaxTypeIndependenceTest",
          definition = function(object, 
              method = c("global", "single-step", "step-down", "discrete"), ...) {

              method <- match.arg(method)
              x <- object@statistic
              RET <- switch(method, 
                   "global" = pvalue(object@distribution, 
                                    object@statistic@teststatistic),
                   "single-step" = singlestep(object, ...),
                   "step-down" = stepdown(object, ...),
                   "discrete" = dbonf(object, ...)
              )
              return(RET)
          }
)

setMethod(f = "pvalue",
          signature = "QuadTypeIndependenceTest",
          definition = function(object, ...) {
              pvalue(object@distribution, object@statistic@teststatistic)
          }
)

### generic method for the permutation distribution from objects
setGeneric("pperm", function(object, q, ...)
    standardGeneric("pperm"))

setMethod(f = "pperm",
          signature = "NullDistribution",
          definition = function(object, q, ...) {
              object@p(q)
          }
)

setMethod(f = "pperm",
          signature = "IndependenceTest",
          definition = function(object, q, ...) {
              pperm(object@distribution, q = q)
          }
)

setMethod(f = "pperm",
          signature = "ScalarIndependenceTest",
          definition = function(object, q, ...) {
              pperm(object@distribution, q)
          }
)

setMethod(f = "pperm",
          signature = "MaxTypeIndependenceTest",
          definition = function(object, q, ...) {
              pperm(object@distribution, q)
          }
)

setMethod(f = "pperm",
          signature = "QuadTypeIndependenceTest",
          definition = function(object, q, ...) {
              pperm(object@distribution, q)
          }
)

### generic method for the permutation distribution from objects
setGeneric("qperm", function(object, p, ...)
    standardGeneric("qperm"))

setMethod(f = "qperm",
          signature = "NullDistribution",
          definition = function(object, p, ...) {
              object@q(p)
          }
)

setMethod(f = "qperm",
          signature = "IndependenceTest",
          definition = function(object, p, ...) {
              qperm(object@distribution, p)
          }
)

setMethod(f = "qperm",
          signature = "ScalarIndependenceTest",
          definition = function(object, p, ...) {
              qperm(object@distribution, p)
          }
)

setMethod(f = "qperm",
          signature = "MaxTypeIndependenceTest",
          definition = function(object, p, ...) {
              qperm(object@distribution, p)
          }
)

setMethod(f = "qperm",
          signature = "QuadTypeIndependenceTest",
          definition = function(object, p, ...) {
              qperm(object@distribution, p)
          }
)

### generic method for the permutation distribution from objects
setGeneric("dperm", function(object, x, ...)
    standardGeneric("dperm"))

setMethod(f = "dperm",
          signature = "NullDistribution",
          definition = function(object, x, ...) {
              object@d(x)
          }
)

setMethod(f = "dperm",
          signature = "IndependenceTest",
          definition = function(object, x, ...) {
              dperm(object@distribution, x)
          }
)

setMethod(f = "dperm",
          signature = "ScalarIndependenceTest",
          definition = function(object, x, ...) {
              dperm(object@distribution, x)
          }
)

setMethod(f = "dperm",
          signature = "MaxTypeIndependenceTest",
          definition = function(object, x, ...) {
              dperm(object@distribution, x)
          }
)

setMethod(f = "dperm",
          signature = "QuadTypeIndependenceTest",
          definition = function(object, x, ...) {
              dperm(object@distribution, x)
          }
)

### generic method for the permutation distribution from objects
setGeneric("support", function(object, ...)
    standardGeneric("support"))

setMethod(f = "support",
          signature = "NullDistribution",
          definition = function(object, ...) {
              object@support(...)
          }
)

setMethod(f = "support",
          signature = "IndependenceTest",
          definition = function(object, ...) {
              support(object@distribution, ...)
          }
)

setMethod(f = "support",
          signature = "ScalarIndependenceTest",
          definition = function(object, ...) {
              support(object@distribution, ...)
          }
)

setMethod(f = "support",
          signature = "MaxTypeIndependenceTest",
          definition = function(object, ...) {
              support(object@distribution, ...)
          }
)

setMethod(f = "support",
          signature = "QuadTypeIndependenceTest",
          definition = function(object, ...) {
              support(object@distribution, ...)
          }
)

### generic method for extracting statistics from objects
setGeneric("statistic", function(object, 
    type = c("test", "linear", "standardized"), ...) 
    standardGeneric("statistic")
)

setMethod(f = "statistic",
          signature = "IndependenceTest",
          definition = function(object, 
              type = c("test", "linear", "standardized"), ...) {
              ### <FIXME> it is not sure that object@statistic exists! </FIXME>
              object <- object@statistic
              type <- match.arg(type)
              nc <- ncol(object@ytrans)
              nr <- ncol(object@xtrans)
              dn <- list(colnames(object@xtrans), 
                         colnames(object@ytrans))
              if (object@has_scores) {
                  if (object@xordinal) {
                      x <- object@x[[1]]
                      nr <- 1
                      dn[[1]] <- paste(attr(x, "scores"), "*", 
                                       abbreviate(levels(x)), 
                                       collapse = " + ", sep = "")
                  }
                  if (object@yordinal) {
                      y <- object@y[[1]]
                      nc <- 1
                      dn[[2]] <- paste(attr(y, "scores"), "*", 
                                       abbreviate(levels(y)), 
                                       collapse = " + ", sep = "")
                  }
              }
              switch(type, "test" = object@teststatistic,
                           "linear" = matrix(object@linearstatistic, 
                                             nrow = nr, ncol = nc, dimnames = dn),
                           "standardized" = matrix(object@standardizedlinearstatistic, 
                                                   nrow = nr, ncol = nc, dimnames = dn)
                           )
      }
)

### generic method for extracting expectations from objects
setGeneric("expectation", function(object, ...) 
    standardGeneric("expectation")
)

setMethod(f = "expectation",
          signature = "IndependenceTest",
          definition = function(object, ...) 
              expectation(object@statistic, ...)
)

setMethod(f = "expectation",
          signature = "IndependenceTestStatistic",
          definition = function(object, ...) 
              object@expectation
)


setMethod(f = "expectation",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, ...) 
              object@expectation
)

setMethod(f = "expectation",
          signature = "MaxTypeIndependenceTestStatistic",
          definition = function(object, ...) {

#              x <- object@expectation
#              matrix(x, nrow = ncol(object@xtrans),
#                        ncol = ncol(object@ytrans),
#                        dimnames = list(colnames(object@xtrans), 
#                                        colnames(object@ytrans)))
               object@expectation
          }
)

setMethod(f = "expectation",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, ...) {

#              x <- object@expectation
#              matrix(x, nrow = ncol(object@xtrans),
#                        ncol = ncol(object@ytrans),
#                        dimnames = list(colnames(object@xtrans), 
#                                        colnames(object@ytrans)))
              object@expectation
          }
)

### generic method for extracting the covariance matrix from objects
setGeneric("covariance", function(object, ...) 
    standardGeneric("covariance")
)

setMethod(f = "covariance",
          signature = "CovarianceMatrix",
          definition = function(object, ...) 
              object@covariance
)

setMethod(f = "covariance",
          signature = "IndependenceTest",
          definition = function(object, ...) 
              covariance(object@statistic, ...)
)

setMethod(f = "covariance",
          signature = "IndependenceTestStatistic",
          definition = function(object, ...) 
              covariance(object@covariance)
)

setMethod(f = "covariance",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, ...) 
              covariance(object@covariance)
)

setMethod(f = "covariance",
          signature = "MaxTypeIndependenceTestStatistic",
          definition = function(object, ...)
              covariance(object@covariance)
)

setMethod(f = "covariance",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, ...)
              covariance(object@covariance)
)

### generic method for extracting the variances 
setGeneric("variance", function(object, ...) 
    standardGeneric("variance")
)

setMethod(f = "variance",
          signature = "Variance",
          definition = function(object, ...) 
              object@variance
)

setMethod(f = "variance",
          signature = "CovarianceMatrix",
          definition = function(object, ...) 
              diag(object@covariance)
)


setMethod(f = "variance",
          signature = "IndependenceTest",
          definition = function(object, ...) 
              variance(object@statistic, ...)
)

setMethod(f = "variance",
          signature = "IndependenceTestStatistic",
          definition = function(object, ...) 
              variance(object@covariance)
)

setMethod(f = "variance",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, ...) 
              variance(object@covariance)
)

setMethod(f = "variance",
          signature = "MaxTypeIndependenceTestStatistic",
          definition = function(object, ...)
              variance(object@covariance)
)

setMethod(f = "variance",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, ...)
              variance(object@covariance)
)

