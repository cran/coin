### <DEPRECATED>
### Note: The "CovarianceMatrix", "Variance" and "VarCovar" classes were
### deprecated in 1.4-0.  To be removed in 2.0-0.
### new("CovarianceMatrix", ...)
setMethod("initialize",
    signature = "CovarianceMatrix",
    definition = function(.Object, covariance, ...) {
        callNextMethod(.Object, covariance = covariance, ...)
    }
)

### new("Variance", ...)
setMethod("initialize",
    signature = "Variance",
    definition = function(.Object, variance, ...) {
        callNextMethod(.Object, variance = variance, ...)
    }
)
### </DEPRECATED>

### new("IndependenceProblem", ...)
### initialized data
setMethod("initialize",
    signature = "IndependenceProblem",
    definition = function(.Object, x, y, block = NULL, weights = NULL, ...) {

        if (NROW(x) == 0L && NROW(y) == 0L)
            stop(sQuote("x"), " and ", sQuote("y"),
                 " do not contain data")
        if (length(x) == 0L) {
            dn <- dimnames(x)
            x <- data.frame(x = rep.int(1L, nrow(x)))
            dimnames(x) <- dn
        }
        if (anyNA(x))
            stop(sQuote("x"), " contains missing values")
        if (anyNA(y))
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && !is.factor(block))
            stop(sQuote("block"), " is not a factor")
        if (!is.null(block) && anyNA(block))
            stop(sQuote("block"), " contains missing values")
        if (!is.null(weights) && anyNA(weights))
            stop(sQuote("weights"), " contains missing values")

        .Object@x <- droplevels(x)
        .Object@y <- droplevels(y)
        .Object@block <- if (is.null(block))
                             factor(rep.int(0L, nrow(x)))
                         else {
                             blockname <- attr(block, "blockname", exact = TRUE)
                             block <- droplevels(block)
                             if (!is.null(blockname))
                                 attr(block, "blockname") <- blockname
                             if (any(table(block) < 2L))
                                 stop(sQuote("block"), " contains levels with",
                                      " less than two observations")
                             block
                         }
        .Object@weights <- if (is.null(weights))
                               rep.int(1.0, nrow(x))
                           else
                               as.double(weights)

        if (!validObject(.Object))
            stop("not a valid object of class ", dQuote("IndependenceProblem"))

        .Object
    }
)

### new("IndependenceTestProblem", ...)
### set up test problem, i.e., transformations of the data
setMethod("initialize",
    signature = "IndependenceTestProblem",
    definition = function(.Object, object, xtrafo = trafo, ytrafo = trafo, ...) {

        if (!inherits(object, "IndependenceProblem"))
            stop(sQuote("object"), " is not of class ",
                 dQuote("IndependenceProblem"))

        tr <- check_trafo(xtrafo(object@x), ytrafo(object@y))

        .Object <- copyslots(object, .Object)
        .Object@xtrans <- tr$xtrafo
        .Object@ytrans <- tr$ytrafo
        .Object@xtrafo <- xtrafo
        .Object@ytrafo <- ytrafo

        .Object
    }
)

### new("IndependenceLinearStatistic", ...)
### compute linear statistics and their expectation / covariance matrix
setMethod("initialize",
    signature = "IndependenceLinearStatistic",
    definition = function(.Object, object, ...) {

        if (!inherits(object, "IndependenceTestProblem"))
            stop(sQuote("object"), " is not of class ",
                 dQuote("IndependenceTestProblem"))

        block <- object@block
        r <- nlevels(block)

        if (r == 1) {
            ecs <- .Call(R_ExpectationCovarianceStatistic,
                         object@xtrans,
                         object@ytrans,
                         object@weights,
                         integer(0), integer(0), 0L, sqrt_eps)
            linearstatistic <- as.matrix(ecs$LinearStatistic)
            expectation <-  as.matrix(ecs$Expectation)
            covariance <- as.matrix(ecs$Covariance)
        } else {
            ytrans <- object@ytrans
            xtrans <- object@xtrans
            weights <- object@weights
            pq <- ncol(xtrans) * ncol(ytrans)

            linearstatistic <- matrix(NA_real_, nrow = pq, ncol = r)
            expectation <- matrix(NA_real_, nrow = pq, ncol = r)
            covariance <- matrix(NA_real_, nrow = pq * (pq + 1) / 2, ncol = r)

            bl <- levels(block)
            for (i in seq_len(r)) {
                block_i <- block == bl[i]
                ecs <- .Call(R_ExpectationCovarianceStatistic,
                             xtrans[block_i,, drop = FALSE],
                             ytrans[block_i,, drop = FALSE],
                             weights[block_i],
                             integer(0), integer(0), 0L, sqrt_eps)
                linearstatistic[, i] <- ecs$LinearStatistic
                expectation[, i] <- ecs$Expectation
                covariance[, i] <- ecs$Covariance
            }
        }

        .Object <- copyslots(object, .Object)
        .Object@linearstatistic <- linearstatistic
        .Object@expectation <- expectation
        .Object@covariance <- covariance

        .Object
    }
)

### compute standardized linear statistics
setMethod("initialize",
    signature = "IndependenceTestStatistic",
    definition = function(.Object, object, ...) {

        if (!inherits(object, "IndependenceLinearStatistic"))
            stop(sQuote("object"), " is not of class ",
                 dQuote("IndependenceLinearStatistic"))

        variance <- .variance(object, partial = FALSE)

        .Object <- copyslots(object, .Object)
        .Object@standardizedlinearstatistic <-
            as.vector(.centeredlinearstatistic(object, partial = FALSE) /
                      sqrt(variance))

        if (any(variance < sqrt_eps))
            warning("The conditional covariance matrix has ",
                    "zero diagonal elements")

        .Object
    }
)

### new("ScalarIndependenceTestStatistic", ...)
### the basis of well known univariate tests
setMethod("initialize",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(.Object, object,
        alternative = c("two.sided", "less", "greater"), paired = FALSE, ...) {

        .Object <- callNextMethod(.Object, object)
        .Object@teststatistic <- .Object@standardizedlinearstatistic
        .Object@alternative <- match.arg(alternative)
        .Object@paired <- paired

        .Object
    }
)

### new("MaxTypeIndependenceTestStatistic", ...)
setMethod("initialize",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(.Object, object,
        alternative = c("two.sided", "less", "greater"), ...) {

        .Object <- callNextMethod(.Object, object)
        .Object@teststatistic <-
            switch(alternative,
                "less"      = min(.Object@standardizedlinearstatistic),
                "greater"   = max(.Object@standardizedlinearstatistic),
                "two.sided" = max(abs(.Object@standardizedlinearstatistic))
            )
        .Object@alternative <- match.arg(alternative)

        .Object
    }
)

### new("QuadTypeIndependenceTestStatistic", ...)
setMethod("initialize",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(.Object, object, paired = FALSE, ...) {

        covarianceplus <- .covariance(object, invert = TRUE, partial = FALSE)

        .Object <- callNextMethod(.Object, object)
        .Object@teststatistic <-
            .Call(R_quadform,
                  .linearstatistic(object, partial = FALSE),
                  .expectation(object, partial = FALSE),
                  covarianceplus)
        .Object@df <- attr(covarianceplus, "rank")
        .Object@covarianceplus <- as.vector(covarianceplus)
        .Object@paired <- paired

        .Object
    }
)

### new("SymmetryProblem", ...)
### initialized data
setMethod("initialize",
    signature = "SymmetryProblem",
    definition = function(.Object, x, y, block = NULL, weights = NULL, ...) {

        if (anyNA(x))
            stop(sQuote("x"), " contains missing values")
        if (!is.factor(x[[1L]]) || length(unique(table(x[[1L]]))) != 1L)
            stop(sQuote("x"), " is not a balanced factor")
        if (anyNA(y))
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && anyNA(y))
            stop(sQuote("block"), " contains missing values")

        .Object@x <- x
        .Object@y <- y
        .Object@block <- if (is.null(block))
                             factor(rep.int(seq_len(nrow(x) / nlevels(x[[1L]])),
                                            nlevels(x[[1L]])))
                         else
                             block
        .Object@weights <- if (is.null(weights))
                               rep.int(1.0, nrow(x))
                           else
                               as.double(weights)

        if (!validObject(.Object))
            stop("not a valid object of class ", sQuote("SymmetryProblem"))

        .Object
    }
)
