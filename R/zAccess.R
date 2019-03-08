### generic method for extracting p-values from objects
setGeneric("pvalue",
    function(object, ...) {
        standardGeneric("pvalue")
    }
)

### <DEPRECATED>
### The "PValue" class was deprecated in 1.3-0 and at the same time this
### method was added as a temporary solution.  To be removed in 2.0-0.
setMethod("pvalue",
    signature = "PValue",
    definition = function(object, q, ...) {
        RET <- object@pvalue(q)
        class(RET) <- "pvalue"
        RET
    }
)
### </DEPRECATED>

setMethod("pvalue",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        RET <- object@pvalue(q)
        class(RET) <- "pvalue"
        RET
    }
)

setMethod("pvalue",
    signature = "ApproxNullDistribution",
    definition = function(object, q, ...) {
        RET <- callNextMethod(object, q, ...)
        attr(RET, "nresample") <- object@nresample
        RET
    }
)

setMethod("pvalue",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, object@statistic@teststatistic, ...)
    }
)

setMethod("pvalue",
    signature = "MaxTypeIndependenceTest",
    definition = function(object,
        method = c("global", "single-step", "step-down", "unadjusted"),
###     combinations = c("free", "restricted"), # placeholder
        distribution = c("joint", "marginal"),
        type = c("Bonferroni", "Sidak"), ...) {
            method <- match.arg(method,
                          choices = c("global", "single-step", "step-down",
                                      "unadjusted", "discrete"),
                          several.ok = TRUE)[1]
            if (method == "discrete")
                stop(sQuote(paste("method =", dQuote(method))),
                        " is defunct; see ", sQuote("?pvalue"))
            distribution <- match.arg(distribution)
            type <- match.arg(type)

            C <- attr(object@statistic@xtrans, "contrast")
            if (!is.null(C) && method != "global")
                warning("p-values may be incorrect due to violation\n",
                        "  of the subset pivotality condition")
            ## NOTE: Two ^^ spaces needed for correct rendering

            if (method == "global")
                callNextMethod(object, ...)
            else if (method == "single-step") {
                if (distribution == "joint")
                    joint(object, stepdown = FALSE, ...)
                else {
                    if (type == "Bonferroni")
                        marginal(object, stepdown = FALSE,
                                 bonferroni = TRUE, ...)
                    else
                        marginal(object, stepdown = FALSE,
                                 bonferroni = FALSE, ...)
                }
            } else if (method == "step-down") {
                if (distribution == "joint")
                    joint(object, stepdown = TRUE, ...)
                else {
                    if (type == "Bonferroni")
                        marginal(object, stepdown = TRUE,
                                 bonferroni = TRUE, ...)
                    else
                        marginal(object, stepdown = TRUE,
                                 bonferroni = FALSE, ...)
                }
            }
            else
                unadjusted(object, ...)
    }
)


### generic method for extracting mid-p-values from objects
setGeneric("midpvalue",
    function(object, ...) {
        standardGeneric("midpvalue")
    }
)

setMethod("midpvalue",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        RET <- object@midpvalue(q)
        class(RET) <- "pvalue"
        RET
    }
)

setMethod("midpvalue",
    signature = "ApproxNullDistribution",
    definition = function(object, q, ...) {
        RET <- callNextMethod(object, q, ...)
        attr(RET, "nresample") <- object@nresample
        RET
    }
)

setMethod("midpvalue",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, object@statistic@teststatistic, ...)
    }
)


### generic method for extracting p-value intervals from objects
setGeneric("pvalue_interval",
    function(object, ...) {
        standardGeneric("pvalue_interval")
    }
)

setMethod("pvalue_interval",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        object@pvalueinterval(q)
    }
)

setMethod("pvalue_interval",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, object@statistic@teststatistic, ...)
    }
)


### generic method for extracting size from objects
setGeneric("size",
    function(object, ...) {
        standardGeneric("size")
    }
)

setMethod("size",
    signature = "NullDistribution",
    definition = function(object,
        alpha, type = c("p-value", "mid-p-value"), ...) {
            type <- match.arg(type)
            object@size(alpha, type)
    }
)

setMethod("size",
    signature = "IndependenceTest",
    definition = function(object,
        alpha, type = c("p-value", "mid-p-value"), ...) {
            callGeneric(object@distribution, alpha, type, ...)
    }
)


### generic method for the permutation distribution from objects
setGeneric("dperm",
    function(object, x, ...) {
        standardGeneric("dperm")
    }
)

setMethod("dperm",
    signature = "NullDistribution",
    definition = function(object, x, ...) {
        object@d(x)
    }
)

setMethod("dperm",
    signature = "IndependenceTest",
    definition = function(object, x, ...) {
        callGeneric(object@distribution, x, ...)
    }
)


### generic method for the permutation distribution from objects
setGeneric("pperm",
    function(object, q, ...) {
        standardGeneric("pperm")
    }
)

setMethod("pperm",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        object@p(q)
    }
)

setMethod("pperm",
    signature = "IndependenceTest",
    definition = function(object, q, ...) {
        callGeneric(object@distribution, q, ...)
    }
)


### generic method for the permutation distribution from objects
setGeneric("qperm",
    function(object, p, ...) {
        standardGeneric("qperm")
    }
)

setMethod("qperm",
    signature = "NullDistribution",
    definition = function(object, p, ...) {
        object@q(p)
    }
)

setMethod("qperm",
    signature = "IndependenceTest",
    definition = function(object, p, ...) {
        callGeneric(object@distribution, p, ...)
    }
)


### generic method for the permutation distribution from objects
setGeneric("rperm",
    function(object, n, ...) {
        standardGeneric("rperm")
    }
)

setMethod("rperm",
    signature = "NullDistribution",
    definition = function(object, n, ...) {
        object@q(runif(n))
    }
)

setMethod("rperm",
    signature = "IndependenceTest",
    definition = function(object, n, ...) {
        callGeneric(object@distribution, n, ...)
    }
)


### generic method for the permutation distribution from objects
setGeneric("support",
    function(object, ...) {
        standardGeneric("support")
    }
)

setMethod("support",
    signature = "NullDistribution",
    definition = function(object, ...) {
        object@support(...)
    }
)

setMethod("support",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, ...)
    }
)


### generic method for extracting statistics from objects
setGeneric("statistic",
    function(object, ...) {
        standardGeneric("statistic")
    }
)

setMethod("statistic",
    signature = "IndependenceLinearStatistic",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"), ...) {
            type <- match.arg(type)
            nr <- ncol(object@xtrans)
            nc <- ncol(object@ytrans)
            dn <- statnames(object)$dimnames
            switch(type,
                "test" = stop(
                    sQuote(paste("type =", dQuote("test"))),
                    " not defined for objects of class ",
                    dQuote("IndependenceLinearStatistic")
                ),
                "linear" = matrix(
                    object@linearstatistic,
                    nrow = nr, ncol = nc, dimnames = dn
                ),
                "centered" = matrix(
                    object@linearstatistic - object@expectation,
                    nrow = nr, ncol = nc, dimnames = dn
                ),
                "standardized" = matrix(
                    object@standardizedlinearstatistic,
                    nrow = nr, ncol = nc, dimnames = dn
                )
            )
    }
)

setMethod("statistic",
    signature = "IndependenceTestStatistic",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"), ...) {
            type <- match.arg(type)
            if (type == "test")
                object@teststatistic
            else
                callNextMethod(object, type, ...)
    }
)

setMethod("statistic",
    signature = "IndependenceTest",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"), ...) {
            callGeneric(object@statistic, type, ...)
    }
)


### generic method for extracting expectations from objects
setGeneric("expectation",
    function(object, ...) {
        standardGeneric("expectation")
    }
)

setMethod("expectation",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        object@expectation
    }
)

setMethod("expectation",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@statistic, ...)
    }
)


### generic method for extracting the covariance matrix from objects
setGeneric("covariance",
    function(object, ...) {
        standardGeneric("covariance")
    }
)

setMethod("covariance",
    signature = "CovarianceMatrix",
    definition = function(object, ...) {
        object@covariance
    }
)

setMethod("covariance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        if (!inherits(object@covariance, "CovarianceMatrix"))
            object <- new("IndependenceLinearStatistic", object, varonly = FALSE)
        callGeneric(object@covariance, ...)
    }
)

setMethod("covariance",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@statistic, ...)
    }
)


### generic method for extracting the variances
setGeneric("variance",
    function(object, ...) {
        standardGeneric("variance")
    }
)

setMethod("variance",
    signature = "Variance",
    definition = function(object, ...) {
        object@variance
    }
)

setMethod("variance",
    signature = "CovarianceMatrix",
    definition = function(object, ...) {
        diag(object@covariance)
    }
)

setMethod("variance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        callGeneric(object@covariance, ...)
    }
)

setMethod("variance",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@statistic, ...)
    }
)
