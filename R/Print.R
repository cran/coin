
varnames <- function(object) {

    x <- object@x
    y <- object@y

    yordered <- sapply(y, is.ordered)
    ynames <- paste(colnames(y), ifelse(yordered, " (ordered)", ""), sep = "", 
                    collapse = ", ")

    if (length(x) == 1) {
        if (is.ordered(x[[1]])) {
            xnames <- paste(colnames(x), " (", paste(levels(x[[1]]), collapse = " < "), 
                            ")", sep = "")
        } else {
            if (is.factor(x[[1]])) {
                xnames <- paste(colnames(x), " (", 
                                paste(levels(x[[1]]), collapse = ", "), ")", sep = "")
            } else {
                xnames <- colnames(x)
            }
        }
    } else {
        xordered <- sapply(x, is.ordered)
        xnames <- paste(colnames(x), ifelse(xordered, "(ordered)", ""), 
                        sep = "", collapse = ", ")
    }

    if (nlevels(object@block) > 1) {
        bn <- attr(object@block, "blockname")
        if (is.null(bn)) bn <- "block"
        xnames <- paste(xnames, paste("\n\t stratified by", bn))
    }

    if (nchar(xnames) > options("width")$width/2) { 
        strg <- paste(ynames, "by\n\t", xnames, collapse = "")
    } else {
        strg <- paste(ynames, "by", xnames, collapse = "")
    }

    return(strg)
}

setMethod(f = "show", signature = "QuadTypeIndependenceTest", 
    definition = function(object) {

        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "chi-squared"
        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact") 

        dataname <- varnames(x@statistic)

        RET <- list(statistic = stat,
                    p.value = dist@pvalue(stat),
                    data.name = dataname,
                    method = paste(distname, x@method))
        if (length(dist@parameters) == 1 && 
                   names(dist@parameters) == "df") {
            RET$parameter <- dist@parameters[[1]]
            names(RET$parameter) <- "df"
        }
        if (length(x@estimates) > 0)
            RET <- c(RET, x@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod(f = "show", signature = "MaxTypeIndependenceTest",
    definition = function(object) {

        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "maxT"
        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        dataname <- varnames(x@statistic)

        RET <- list(statistic = stat,
                    p.value = x@distribution@pvalue(stat),
                    data.name = dataname,
                    method = paste(distname, x@method))
        if (length(x@estimates) > 0)
            RET <- c(RET, x@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod(f = "show", signature = "ScalarIndependenceTest",
    definition = function(object) {
          
        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "Z"

        dataname <- varnames(x@statistic)

        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(statistic = stat,
                    p.value = x@distribution@pvalue(stat),
                    alternative = x@statistic@alternative,
                    data.name = dataname,
                    method = paste(distname, x@method))
        if (length(x@nullvalue > 0)) 
            RET$null.value = c(mu = x@nullvalue)
        if (length(x@estimates) > 0)
            RET <- c(RET, x@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod(f = "show", signature = "IndependenceTest",
    definition = function(object) {
          
        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "c"

        dataname <- varnames(x@statistic)

        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(statistic = stat,
                    p.value = x@distribution@pvalue(stat),
                    data.name = dataname,
                    method = paste(distname, x@method))
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod(f = "show", signature = "ScalarIndependenceTestConfint",
    definition = function(object) {

        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "Z"
        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        dataname <- varnames(x@statistic)

        ci <- confint(object, level = object@conf.level)

        RET <- list(statistic = stat,
                    p.value = x@distribution@pvalue(stat),
                    alternative = x@statistic@alternative,
                    method = paste(distname, x@method),
                    data.name = dataname,
                    conf.int = ci$conf.int,
                    estimate = ci$estimate)
        if (length(x@nullvalue))
            RET$null.value = c(mu = x@nullvalue)
        if (length(x@estimates))
            RET <- c(RET, x@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

print.ci <- function(x, ...) {

    if(!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")),
            "percent confidence interval:\n",
            format(c(x$conf.int[1], x$conf.int[2])), "\n")
    }
    if(!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}

print.MCp <- function(x, ...) {

    p <- x
    attributes(p) <- NULL
    print(p)
    ci <- list(conf.int = attr(x, "conf.int"))
    class(ci) <- "ci"
    print(ci)
}
