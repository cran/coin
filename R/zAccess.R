
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
              pvalue(object@distribution, object@statistic@teststatistic)
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

### generic method for extracting statistics from objects
setGeneric("statistic", function(object, 
    type = c("test", "linear", "standardized"), ...) 
    standardGeneric("statistic")
)

setMethod(f = "statistic",
          signature = "IndependenceTest",
          definition = function(object, 
              type = c("test", "linear", "standardized"), ...)
              statistic(object@statistic, type = type)
)

setMethod(f = "statistic",
          signature = "IndependenceTestStatistic",
          definition = function(object, 
              type = c("test", "linear", "standardized"), ...) {
              nc <- ncol(object@ytrans)
              nr <- ncol(object@xtrans)
              type <- match.arg(type)
              dn <- statnames(object)$dimnames
              switch(type, "test" = object@teststatistic,
                           "linear" = matrix(object@linearstatistic, 
                                             nrow = nr, ncol = nc, dimnames = dn),
                           "standardized" = matrix(object@standardizedlinearstatistic, 
                                                   nrow = nr, ncol = nc, dimnames = dn)
                           )
      }
)

setMethod(f = "statistic",
          signature = "IndependenceLinearStatistic",
          definition = function(object, 
              type = c("test", "linear", "standardized"), ...) {
              nc <- ncol(object@ytrans)
              nr <- ncol(object@xtrans)
              type <- match.arg(type)
              dn <- statnames(object)$dimnames
              switch(type, "test" = stop("type = test not defined for IndependenceLinearStatistic"),
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
          signature = "IndependenceLinearStatistic",
          definition = function(object, ...) 
              object@expectation
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
          signature = "IndependenceLinearStatistic",
          definition = function(object, ...) {
              if (!extends(class(object@covariance), "CovarianceMatrix"))
                  return(covariance(new("IndependenceTestStatistic", 
                                        object, varonly = FALSE)))
              covariance(object@covariance)
          }
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
          signature = "IndependenceLinearStatistic",
          definition = function(object, ...) 
              variance(object@covariance)
)
