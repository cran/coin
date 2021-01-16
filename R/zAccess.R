### methods for extracting p-values
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
        class(RET) <- c("pvalue", "numeric")
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


### methods for extracting mid-p-values
setGeneric("midpvalue",
    function(object, ...) {
        standardGeneric("midpvalue")
    }
)

setMethod("midpvalue",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        RET <- object@midpvalue(q)
        class(RET) <- c("pvalue", "numeric")
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


### methods for extracting p-value intervals
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


### methods for extracting test size
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


### methods for extracting the density function
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


### methods for extracting the distribution function
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


### methods for extracting the quantile function
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


### methods for extracting random deviates
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


### methods for extracting the support
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


### *internal* methods for extracting the linear statistic
setGeneric(".linearstatistic",
    function(object, ...) {
        standardGeneric(".linearstatistic")
    }
)

setMethod(".linearstatistic",
    signature = "IndependenceLinearStatistic",
    definition = function(object, partial, ...) {
        RET <- object@linearstatistic
        r <- ncol(RET)
        if (r > 1 && !partial)
            RET <- as.matrix(rowSums(RET))
        RET
    }
)


### *internal* methods for extracting the centered linear statistic
setGeneric(".centeredlinearstatistic",
    function(object, ...) {
        standardGeneric(".centeredlinearstatistic")
    }
)

setMethod(".centeredlinearstatistic",
    signature = "IndependenceLinearStatistic",
    definition = function(object, partial, ...) {
        RET <- object@linearstatistic - object@expectation
        r <- ncol(RET)
        if (r > 1 && !partial)
            RET <- as.matrix(rowSums(RET))
        RET
    }
)


### methods for extracting test statistics etc.
setGeneric("statistic",
    function(object, ...) {
        standardGeneric("statistic")
    }
)

setMethod("statistic",
    signature = "IndependenceLinearStatistic",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"),
        partial = FALSE, ...) {
            type <- match.arg(type)
            p <- ncol(object@xtrans)
            q <- ncol(object@ytrans)
            r <- nlevels(object@block)
            dn <- statnames(object)$dimnames
            RET <- switch(type,
                       "test" = stop(
                           sQuote(paste("type =", dQuote("test"))),
                           " not defined for objects of class ",
                           dQuote("IndependenceLinearStatistic")
                       ),
                       "linear" = {
                           .linearstatistic(object, partial, ...)
                       },
                       "centered" = {
                           .centeredlinearstatistic(object, partial, ...)
                       },
                       "standardized" = {
                           .centeredlinearstatistic(object, partial, ...) /
                               sqrt(.variance(object, partial, ...))
                       }
                   )
            if (r > 1 && partial)
                setAttributes(RET, list(dim = c(p, q, r), dimnames = dn))
            else
                setAttributes(RET, list(dim = c(p, q), dimnames = dn[-3]))
    }
)

setMethod("statistic",
    signature = "IndependenceTestStatistic",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"),
        partial = FALSE, ...) {
            type <- match.arg(type)
            if (type == "test")
                object@teststatistic
            else if (type == "standardized" && !partial) {
                    p <- ncol(object@xtrans)
                    q <- ncol(object@ytrans)
                    dn <- statnames(object)$dimnames
                    matrix(
                        object@standardizedlinearstatistic,
                        nrow = p, ncol = q, dimnames = dn[-3]
                    )
            } else
                callNextMethod(object, type, partial, ...)
    }
)

setMethod("statistic",
    signature = "IndependenceTest",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"),
        partial = FALSE, ...) {
            callGeneric(object@statistic, type, partial, ...)
    }
)


### *internal* methods for extracting expectations
setGeneric(".expectation",
    function(object, ...) {
        standardGeneric(".expectation")
    }
)

setMethod(".expectation",
    signature = "IndependenceLinearStatistic",
    definition = function(object, partial, ...) {
        RET <- object@expectation
        r <- ncol(RET)
        if (r > 1 && !partial)
            RET <- as.matrix(rowSums(RET))
        RET
    }
)


### methods for extracting expectations
setGeneric("expectation",
    function(object, ...) {
        standardGeneric("expectation")
    }
)

setMethod("expectation",
    signature = "IndependenceLinearStatistic",
    definition = function(object, partial = FALSE, ...) {
        RET <- .expectation(object, partial, ...)
        p <- ncol(object@xtrans)
        q <- ncol(object@ytrans)
        r <- ncol(RET)
        dn <- statnames(object)$dimnames
        if (r > 1 && partial)
            setAttributes(RET, list(dim = c(p, q, r), dimnames = dn))
        else
            setAttributes(RET, list(dim = c(p, q), dimnames = dn[-3]))
    }
)

setMethod("expectation",
    signature = "IndependenceTest",
    definition = function(object, partial = FALSE, ...) {
        callGeneric(object@statistic, partial, ...)
    }
)


### *internal* methods for extracting the covariances
setGeneric(".covariance",
    function(object, ...) {
        standardGeneric(".covariance")
    }
)

setMethod(".covariance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, invert, partial, ...) {
        RET <- object@covariance
        r <- ncol(RET)
        if (r > 1 && partial) {
            if (invert) {
                mp_rank <- integer(r)
                for (i in seq_len(r)) {
                    mp <- .Call(R_MPinv_sym, RET[, i], 0L, sqrt_eps)
                    RET[, i] <- mp$MPinv
                    mp_rank[i] <- mp$rank
                }
            }
            RET
        } else {
            if (r > 1)
                RET <- rowSums(RET)
            if (invert) {
                mp <- .Call(R_MPinv_sym, RET, 0L, sqrt_eps)
                RET <- mp$MPinv
                mp_rank <- mp$rank
            }
            RET <- as.matrix(RET)
        }

        if (invert)
            attr(RET, "rank") <- mp_rank
        RET
    }
)


### methods for extracting the covariance matrix
setGeneric("covariance",
    function(object, ...) {
        standardGeneric("covariance")
    }
)

### <DEPRECATED>
### Note: The "CovarianceMatrix", "Variance" and "VarCovar" classes were
### deprecated in 1.4-0.  To be removed in 2.0-0.
setMethod("covariance",
    signature = "CovarianceMatrix",
    definition = function(object, ...) {
        object@covariance
    }
)
### </DEPRECATED>

setMethod("covariance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, invert = FALSE, partial = FALSE, ...) {
        RET <- .covariance(object, invert, partial, ...)
        pq <- nrow(object@linearstatistic)
        r <- ncol(RET)
        nm <- statnames(object)
        RET <- unlist(lapply(seq_len(r), function(i)
            .Call(R_unpack_sym, RET[, i], NULL, 0L)))
        if (r > 1 && partial)
            setAttributes(RET, list(dim = c(pq, pq, r),
                                    dimnames = list(nm$names, nm$names,
                                                    nm$dimnames[[3]])))
        else
            setAttributes(RET, list(dim = c(pq, pq),
                                    dimnames = list(nm$names, nm$names)))
    }
)

setMethod("covariance",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, invert = FALSE, partial = FALSE, ...) {
        nm <- statnames(object)$names
        if (invert && !partial) {
            .Call(R_unpack_sym, object@covarianceplus, nm, 0L)
        } else
            callNextMethod(object, invert, partial, ...)
    }
)

setMethod("covariance",
    signature = "IndependenceTest",
    definition = function(object, invert = FALSE, partial = FALSE, ...) {
        callGeneric(object@statistic, invert, partial, ...)
    }
)


### *internal* methods for extracting the variances
setGeneric(".variance",
    function(object, ...) {
        standardGeneric(".variance")
    }
)

setMethod(".variance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, partial, ...) {
        RET <- object@covariance
        r <- ncol(RET)
        RET <- do.call("cbind", lapply(seq_len(r), function(i)
            .Call(R_unpack_sym, RET[, i], NULL, 1L)))
        if (r > 1 && !partial)
            RET <- as.matrix(rowSums(RET))
        RET
    }
)


### methods for extracting the variances
setGeneric("variance",
    function(object, ...) {
        standardGeneric("variance")
    }
)

### <DEPRECATED>
### Note: The "CovarianceMatrix", "Variance" and "VarCovar" classes were
### deprecated in 1.4-0.  To be removed in 2.0-0.
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
### </DEPRECATED>

setMethod("variance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, partial = FALSE, ...) {
        RET <- .variance(object, partial, ...)
        p <- ncol(object@xtrans)
        q <- ncol(object@ytrans)
        r <- ncol(RET)
        dn <- statnames(object)$dimnames
        if (r > 1 && partial)
            setAttributes(RET, list(dim = c(p, q, r), dimnames = dn))
        else
            setAttributes(RET, list(dim = c(p, q), dimnames = dn[-3]))
    }
)

setMethod("variance",
    signature = "IndependenceTest",
    definition = function(object, partial = FALSE, ...) {
        callGeneric(object@statistic, partial, ...)
    }
)
