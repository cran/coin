
### new("IndependenceProblem", ...)
### initialized data
setMethod(f = "initialize", 
    signature = "IndependenceProblem", 
    definition = function(.Object, x, y, block = NULL, weights = NULL) {

        if (any(is.na(x))) 
            stop(sQuote("x"), " contains missing values")
        if (any(is.na(y))) 
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && any(is.na(y))) 
            stop(sQuote("block"), " contains missing values")
        .Object@x <- x
        .Object@y <- y
        if (is.null(block)) {
            .Object@block <- factor(rep(0, nrow(x)))
        } else {
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
        xfact <- sapply(x, is.factor)
        yfact <- sapply(y, is.factor)

        tr <- check_trafo(xtrafo(x), ytrafo(y))
        .Object@xtrans <- tr$xtrafo
        .Object@ytrans <- tr$ytrafo
        .Object@xtrafo <- xtrafo
        .Object@ytrafo <- ytrafo
        p <- ncol(.Object@xtrans)
        q <- ncol(.Object@ytrans)
        .Object@scores <- diag(p * q)

        if (((ncol(x) > 1 && ncol(tr$xtrafo) > 1) || 
             (ncol(y) > 1 && ncol(tr$ytrafo) > 1)) && 
             any(xfact || yfact)) {
            colnames(.Object@xtrans) <- paste(
                rep(colnames(x), table(attr(.Object@xtrans, "assign"))), 
                    colnames(.Object@xtrans), sep = ".")
            colnames(.Object@xtrans)[attr(.Object@xtrans, "assign") 
                %in% which(!xfact)] <- colnames(x)[!xfact]
            colnames(.Object@ytrans) <- paste(
                rep(colnames(y), table(attr(.Object@ytrans, "assign"))), 
                    colnames(.Object@ytrans), sep = ".")
            colnames(.Object@ytrans)[attr(.Object@ytrans, "assign") 
                %in% which(!yfact)] <- colnames(y)[!yfact]
        }

        ### check if scores are attached
        ### <FIXME> more careful checks!
        xORDINAL <- (ncol(x) == 1 && is.ordered(x[[1]])) && 
                    (nlevels(x[[1]]) > 2)
        yORDINAL <- (ncol(y) == 1 && is.ordered(y[[1]])) && 
                    (nlevels(y[[1]]) > 2) 
        .Object@has_scores <- xORDINAL || yORDINAL

        if (xORDINAL) {
            if (is.null(attr(x[[1]], "scores")) && is.null(xscores))
                xscores <- 1:nlevels(x[[1]])
            else {
                if (is.null(xscores)) 
                    xscores <- attr(x[[1]], "scores")
            }
            if (length(xscores) != ncol(.Object@xtrans)) 
                stop(sQuote("xscores"), " don't match")
        }
        if (yORDINAL) {
            if (is.null(attr(y[[1]], "scores")))
                yscores <-  1:nlevels(y[[1]])
            else {
                if (is.null(yscores))
                    yscores <- attr(y[[1]], "scores")
            }
            if (length(yscores) != ncol(.Object@ytrans)) 
                stop(sQuote("yscores"), " don't match")
        }
        if (xORDINAL && !yORDINAL)
            .Object@scores <- .Call("R_scmatleft",
                as.double(xscores), as.integer(p * q), 
                PACKAGE = "coin")
        if (xORDINAL && yORDINAL)
            ### grrr: class(kronecker(1:4, 1:3)) == "array"  
            .Object@scores <- matrix(kronecker(yscores, xscores), nrow = 1)
        if (!xORDINAL && yORDINAL)
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
    definition = function(.Object, itp) {

        if (!extends(class(itp), "IndependenceTestProblem"))
            stop("Argument ", sQuote("itp"), " is not of class ",
                  sQuote("IndependenceTestProblem"))

        .Object <- copyslots(itp, .Object)

        xtrans <- itp@xtrans
        ytrans <- itp@ytrans
        weights <- itp@weights
        S <- itp@scores

        .Object@linearstatistic <- drop(S %*% LinearStatistic(xtrans, 
                                                              ytrans, weights))
        
        ### <REMAINDER>
        ### for teststat = "maxtype" and distribution = "approx"
        ### we don't need to covariance matrix but the variances only
        ### </REMAINDER>

        ### possibly stratified by block
        if (nlevels(itp@block) == 1) {
            expcov <- ExpectCovarLinearStatistic(xtrans, ytrans, weights)
            exp <- expcov@expectation
            cov <- expcov@covariance
        } else {
            exp <- 0
            cov <- 0
            for (lev in levels(itp@block)) {
                indx <- (itp@block == lev)
                ec <- ExpectCovarLinearStatistic(xtrans[indx,,drop = FALSE], 
                                                 ytrans[indx,,drop = FALSE], 
                                                 weights[indx])
                exp <- exp + ec@expectation
                cov <- cov + ec@covariance
            }
        }
     
        ### multiply with score matrix
        .Object@expectation <- drop(S %*% exp)
        .Object@covariance <- S %*% cov %*% t(S)
        if (any(diag(.Object@covariance) < sqrt(.Machine$double.eps)))
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

        standstat <- (its@linearstatistic - its@expectation) / 
                     sqrt(its@covariance)
        .Object@teststatistic <- drop(standstat)

        .Object
    }
)

### new("MaxTypeIndependenceTestStatistic", ...)
setMethod(f = "initialize", 
    signature = "MaxTypeIndependenceTestStatistic", 
    definition = function(.Object, its) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)

        standstat <- (its@linearstatistic - its@expectation) / 
                      sqrt(diag(its@covariance))
        .Object@teststatistic <- drop(max(abs(standstat)))
        .Object@standardizedlinearstatistic <- standstat

        .Object
    }
)

### new("QuadTypeIndependenceTestStatistic", ...)
setMethod(f = "initialize", 
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(.Object, its, tol = sqrt(.Machine$double.eps)) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)

        covm <- its@covariance
        mp <- MPinv(covm)
        .Object@covarianceplus <- mp$MPinv
        .Object@df <- mp$rank

        stand <- (its@linearstatistic - its@expectation)
        .Object@teststatistic <- 
            drop(stand %*% .Object@covarianceplus %*% stand)

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
