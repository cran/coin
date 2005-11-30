
### new("CovarianceMatrix", ...)
setMethod(f = "initialize", 
    signature = "CovarianceMatrix",
    definition = function(.Object, x) {
        .Object@covariance <- x
        .Object
    }
)

### new("Variance", ...)
setMethod(f = "initialize", 
    signature = "Variance",
    definition = function(.Object, x) {
        .Object@variance <- x
        .Object
    }
)

### new("IndependenceProblem", ...)
### initialized data
setMethod(f = "initialize", 
    signature = "IndependenceProblem", 
    definition = function(.Object, x, y, block = NULL, weights = NULL) {

        if (any(is.na(x))) 
            stop(sQuote("x"), " contains missing values")
        if (any(is.na(y))) 
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && !is.factor(block))
            stop(sQuote("block"), " is not a factor")
        if (!is.null(block) && any(is.na(block))) 
            stop(sQuote("block"), " contains missing values")
        if (!is.null(weights) && any(is.na(weights))) 
            stop(sQuote("weights"), " contains missing values")
        .Object@x <- x
        .Object@y <- y
        if (is.null(block)) {
            .Object@block <- factor(rep(0, nrow(x)))
        } else {
            if (any(table(block) < 2))
                stop(sQuote("block"), 
                     " contains levels with less than two observations")
            .Object@block <- block
        }
        if (is.null(weights)) {
            .Object@weights <- rep(1.0, nrow(x))
        } else {
            .Object@weights <- as.double(weights)
        }
        if (!validObject(.Object))
            stop("not a valid object of class ",
                 sQuote("IndependenceProblem"))
        .Object
    }
)

### new("IndependenceTestProblem", ...)
### set up test problem, i.e., transformations of the data
setMethod(f = "initialize", 
    signature = "IndependenceTestProblem",
    definition = function(.Object, ip, xtrafo = trafo, ytrafo = trafo, 
                          xscores = NULL, yscores = NULL, ...) {

        if (!extends(class(ip), "IndependenceProblem"))
            stop("Argument ", sQuote("ip"), " is not of class ", 
                  sQuote("IndependenceProblem"))

        .Object <- copyslots(ip, .Object)

        x <- ip@x
        y <- ip@y

        tr <- check_trafo(xtrafo(x), ytrafo(y))
        .Object@xtrans <- tr$xtrafo
        .Object@ytrans <- tr$ytrafo
        .Object@xtrafo <- xtrafo
        .Object@ytrafo <- ytrafo
        p <- ncol(.Object@xtrans)
        q <- ncol(.Object@ytrans)
        .Object@scores <- diag(p * q)

        xORDINAL <- sapply(x, is.ordered)
        yORDINAL <- sapply(y, is.ordered)

        ### <FIXME> implement handling of multiple ordered factors
        if ((any(xORDINAL) && length(xORDINAL) > 1) ||
            (any(yORDINAL) && length(yORDINAL) > 1))
            stop("handling of multiple ordered factors currently not implemented")

        .Object@has_scores <- xORDINAL || yORDINAL
        .Object@xordinal <- any(xORDINAL)
        .Object@yordinal <- any(yORDINAL)

        xscores <- c()
        for (i in 1:ncol(x)) {
            if (is.ordered(x[[i]])) {
                sc <- attr(x[[i]], "scores")
                if (is.null(sc)) sc <- 1:nlevels(x[[i]])
            } else {
                sc <- rep(0, sum(attr(tr$xtrafo, "assign") == i))
            }
            xscores <- c(xscores, sc)
        }

        yscores <- c()
        for (i in 1:ncol(y)) {
            if (is.ordered(y[[i]])) {
                sc <- attr(y[[i]], "scores")
                if (is.null(sc)) sc <- 1:nlevels(y[[i]])
            } else {
                sc <- rep(0, sum(attr(tr$ytrafo, "assign") == i))
            }
            yscores <- c(yscores, sc)
        }

        if (any(xORDINAL) && !any(yORDINAL))
            .Object@scores <- .Call("R_scmatleft",
                as.double(xscores), as.integer(p * q), 
                PACKAGE = "coin")
        if (any(xORDINAL) && any(yORDINAL))
            ### grrr: class(kronecker(1:4, 1:3)) == "array"  
            .Object@scores <- matrix(kronecker(yscores, xscores), nrow = 1)
        if (!any(xORDINAL) && any(yORDINAL))
            .Object@scores <- .Call("R_scmatright", 
                    as.double(yscores), as.integer(p * q),
                    PACKAGE = "coin")
        ### </FIXME>
        .Object
    }
)

### new("IndependenceTestStatistic", ...)
### compute test statistics and their expectation / covariance matrix
setMethod(f = "initialize", 
    signature = "IndependenceTestStatistic", 
    definition = function(.Object, itp, varonly = FALSE) {

        if (!extends(class(itp), "IndependenceTestProblem"))
            stop("Argument ", sQuote("itp"), " is not of class ",
                  sQuote("IndependenceTestProblem"))

        .Object <- copyslots(itp, .Object)

        xtrans <- itp@xtrans
        ytrans <- itp@ytrans
        weights <- itp@weights
        SCORES <- itp@has_scores
        S <- itp@scores
        varonly <- varonly && (!SCORES)

        if (SCORES) {
            .Object@linearstatistic <- drop(S %*% LinearStatistic(xtrans, 
                                       ytrans, weights))
        } else {
            .Object@linearstatistic <- drop(LinearStatistic(xtrans,
                                       ytrans, weights))
        }
        
        ### <REMINDER>
        ### for teststat = "maxtype" and distribution = "approx"
        ### we don't need to covariance matrix but the variances only
        ### </REMINDER>

        ### possibly stratified by block
        if (nlevels(itp@block) == 1) {
            expcov <- ExpectCovarLinearStatistic(xtrans, ytrans, weights,
                                                 varonly = varonly)
            exp <- expcov@expectation
            cov <- expcov@covariance
        } else {
            exp <- 0
            cov <- 0
            for (lev in levels(itp@block)) {
                indx <- (itp@block == lev)
                ec <- ExpectCovarLinearStatistic(xtrans[indx,,drop = FALSE], 
                                                 ytrans[indx,,drop = FALSE], 
                                                 weights[indx],
                                                 varonly = varonly)
                exp <- exp + ec@expectation
                cov <- cov + ec@covariance
            }
        }


        ### multiply with score matrix if necessary
        if (SCORES) {
            .Object@expectation <- drop(S %*% exp)
            .Object@covariance <- new("CovarianceMatrix", S %*% cov %*% t(S))
        } else {
            .Object@expectation <- drop(exp)
            if (varonly) {
                .Object@covariance <- new("Variance", drop(cov))
            } else {
                .Object@covariance <- new("CovarianceMatrix", cov)
            }
        }

        ### pretty names
        nm <- statnames(itp)$names
        names(.Object@expectation) <- nm

        if (extends(class(.Object@covariance), "CovarianceMatrix")) {
                colnames(.Object@covariance@covariance) <- nm
                rownames(.Object@covariance@covariance) <- nm
        }
        if (extends(class(.Object@covariance), "Variance"))
                names(.Object@covariance@variance) <- nm

        if (any(variance(.Object) < eps()))
            warning("The conditional covariance matrix has ",
                    "zero diagonal elements")
        .Object
    }
)

### new("ScalarIndependenceTestStatistic", ...)
### the basis of well known univariate tests
setMethod(f = "initialize", 
    signature = "ScalarIndependenceTestStatistic", 
    definition = function(.Object, its, 
        alternative = c("two.sided", "less", "greater")) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)
        .Object@alternative <- match.arg(alternative)

        standstat <- (its@linearstatistic - expectation(its)) / 
                     sqrt(variance(its))
        .Object@teststatistic <- drop(standstat)
        .Object@standardizedlinearstatistic <- drop(standstat)

        .Object
    }
)

### new("MaxTypeIndependenceTestStatistic", ...)
setMethod(f = "initialize", 
    signature = "MaxTypeIndependenceTestStatistic", 
    definition = function(.Object, its, 
        alternative = c("two.sided", "less", "greater")) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)

        .Object@alternative <- match.arg(alternative)
        standstat <- (its@linearstatistic - expectation(its)) / 
                      sqrt(variance(its))
        .Object@teststatistic <- switch(alternative,
            "less" = drop(min(standstat)),
            "greater" = drop(max(standstat)),
            "two.sided" = drop(max(abs(standstat)))
         )
        .Object@standardizedlinearstatistic <- standstat

        .Object
    }
)

### new("QuadTypeIndependenceTestStatistic", ...)
setMethod(f = "initialize", 
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(.Object, its, ...) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)

        covm <- covariance(its)
        mp <- MPinv(covm, ...)
        .Object@covarianceplus <- mp$MPinv
        .Object@df <- mp$rank

        stand <- (its@linearstatistic - expectation(its))
        .Object@teststatistic <- 
            drop(stand %*% .Object@covarianceplus %*% stand)
        standstat <- (its@linearstatistic - expectation(its)) /
                      sqrt(variance(its))
        .Object@standardizedlinearstatistic <- standstat

        .Object
    }
)

### new("SymmetryProblem", ...)
### initialized data
setMethod(f = "initialize", 
    signature = "SymmetryProblem", 
    definition = function(.Object, x, y, block = NULL, weights = NULL) {

        if (any(is.na(x))) 
            stop(sQuote("x"), " contains missing values")
        if (!is.factor(x[[1]]) || length(unique(table(x[[1]]))) != 1)
            stop(sQuote("x"), " is not a balanced factor")
        if (any(is.na(y))) 
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && any(is.na(y))) 
            stop(sQuote("block"), " contains missing values")
        .Object@x <- x
        .Object@y <- y
        if (is.null(block)) {
            nbl <- nrow(x)/nlevels(x[[1]])
            lindx  <- tapply(1:nrow(x), x[[1]], function(x) x)
            bl <- rep(0, nrow(x))
            for (l in lindx)
                bl[l] <- 1:nbl
            .Object@block <- factor(unlist(bl))
        } else {
            .Object@block <- block
        }
        if (is.null(weights)) {
            .Object@weights <- rep(1.0, nrow(x))
        } else {
            .Object@weights <- as.double(weights)
        }
        .Object
    }
)
