
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
          definition = function(object, adjusted = FALSE, ...) {
              x <- object@statistic
              if (adjusted) {
                  padj <- 1 - sapply(abs(x@standardizedlinearstatistic), pperm, 
                                     object = object, ...)
                  RET <- matrix(padj, nrow = ncol(x@xtrans), 
                                ncol = ncol(x@ytrans), 
                                dimnames = list(colnames(x@xtrans), 
                                                colnames(x@ytrans)))
              } else {
                  RET <- pvalue(object@distribution, 
                                object@statistic@teststatistic)
              }
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
setGeneric("support", function(object, p, ...)
    standardGeneric("support"))

setMethod(f = "support",
          signature = "NullDistribution",
          definition = function(object, p, ...) {
              object@support(p)
          }
)

setMethod(f = "support",
          signature = "IndependenceTest",
          definition = function(object, p, ...) {
              support(object@distribution, p)
          }
)

setMethod(f = "support",
          signature = "ScalarIndependenceTest",
          definition = function(object, p, ...) {
              support(object@distribution, p)
          }
)

setMethod(f = "support",
          signature = "MaxTypeIndependenceTest",
          definition = function(object, p, ...) {
              support(object@distribution, p)
          }
)

setMethod(f = "support",
          signature = "QuadTypeIndependenceTest",
          definition = function(object, p, ...) {
              support(object@distribution, p)
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
              type <- match.arg(type)
              ### <FIXME> it is not sure that object@statistic exists! </FIXME>
              statistic(object@statistic, type = type, ...)
          }
)

setMethod(f = "statistic",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, 
              type = c("test", "linear", "standardized"), ...) {
              type <- match.arg(type)
              if (type == "standardized") type <- "test"
              switch(type, "test" = object@teststatistic,
                           "linear" = {
                               x <- object@linearstatistic
                               matrix(x, nrow = ncol(object@xtrans), 
                                      ncol = ncol(object@ytrans),
                                      dimnames = list(colnames(object@xtrans), 
                                                      colnames(object@ytrans)))
                           })
          }
)

setMethod(f = "statistic",
          signature = "MaxTypeIndependenceTestStatistic",
          definition = function(object, 
              type = c("test", "linear", "standardized"), ...) {
              type <- match.arg(type)
              switch(type, "test" = object@teststatistic,
                           "linear" = {
                               x <- object@linearstatistic
                               matrix(x, nrow = ncol(object@xtrans), 
                                      ncol = ncol(object@ytrans),
                                      dimnames = list(colnames(object@xtrans), 
                                                      colnames(object@ytrans)))
                           },
                           "standardized" = {
                               x <- object@standardizedlinearstatistic
                               matrix(x, nrow = ncol(object@xtrans), 
                                      ncol = ncol(object@ytrans),
                                      dimnames = list(colnames(object@xtrans), 
                                                      colnames(object@ytrans)))
                           })
              }
)

setMethod(f = "statistic",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, type = c("test", "linear", "standardized"), 
              ...) {
              type <- match.arg(type)
              switch(type, "test" = object@teststatistic,
                           "linear" = {
                               x <- object@linearstatistic
                               matrix(x, nrow = ncol(object@xtrans),
                                      ncol = ncol(object@ytrans),
                                      dimnames = list(colnames(object@xtrans), 
                                                      colnames(object@ytrans)))
                           },
                           "standardized" = {
                               matrix((object@linearstatistic - expectation(object)) /
                                     sqrt(variance(object)), 
                                     nrow = ncol(object@xtrans),
                                     ncol = ncol(object@ytrans),
                                     dimnames = list(colnames(object@xtrans),
                                                     colnames(object@ytrans)))
                           }
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

