setMethod("show",
    signature = "IndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "c"),
            p.value = pvalue(object),
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        class(RET) <- "htest2"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "MaxTypeIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "maxT"),
            p.value = pvalue(object),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest2"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "QuadTypeIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "chi-squared"),
            p.value = pvalue(object),
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        parameters <- object@distribution@parameters
        if (length(parameters) == 1 && names(parameters) == "df")
            RET$parameter <- setNames(parameters[[1]], nm = "df")
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest2"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "Z"),
            p.value = pvalue(object),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        nullvalue <- object@nullvalue
        if (length(nullvalue) > 0)
            RET$null.value <- setNames(nullvalue, nm = object@parameter)
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest2"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTestConfint",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )
        ci <- confint(object)

        RET <- list(
            statistic = setNames(statistic(object), nm = "Z"),
            p.value = pvalue(object),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method),
            conf.int = ci$conf.int,
            estimate = ci$estimate
        )
        nullvalue <- object@nullvalue
        if (length(nullvalue) > 0)
            RET$null.value <- setNames(nullvalue, nm = object@parameter)
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest2"
        print(RET)
        invisible(RET)
    }
)

print.htest2 <-
    function(x, digits = getOption("digits"), ...)
{
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")
    out <- character()
    if (!is.null(x$statistic))
        out <- c(out, paste(names(x$statistic), "=",
                            format(x$statistic, digits = max(1L, digits - 2L))))
    if (!is.null(x$parameter))
        out <- c(out, paste(names(x$parameter), "=",
                            format(x$parameter, digits = max(1L, digits - 2L))))
    if (!is.null(x$p.value)) {
        tol <- if (is.null(nresample <- attr(x$p.value, "nresample"))) eps
               else 1 / nresample
        fp <- format.pval(x$p.value, digits = max(1L, digits - 3L), eps = tol)
        out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp
                                       else paste("=", fp)))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1L) {
                alt.char <- switch(x$alternative,
                                two.sided = "not equal to",
                                less      = "less than",
                                greater   = "greater than"
                            )
                cat("true ", names(x$null.value), " is ", alt.char, " ",
                    x$null.value, "\n", sep = "")
            }
            else {
                cat(x$alternative, "\nnull values:\n", sep = "")
                print(x$null.value, digits = digits, ...)
            }
        }
        else cat(x$alternative, "\n", sep = "")
    }
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")),
            " percent confidence interval:\n", " ",
            paste(format(x$conf.int[1:2], digits = digits), collapse = " "),
            "\n", sep = "")
    }
    if (!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, digits = digits, ...)
    }
    cat("\n")
    invisible(x)
}

print.ci <- function(x, ...) {
    if (hasName(x, "conf.int")) {
        cat(format(100 * attr(x$conf.int, "conf.level")),
            "percent confidence interval:\n",
            format(c(x$conf.int[1], x$conf.int[2]), ...), "\n")
    }
    if (hasName(x, "estimate")) {
        cat("sample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}

format.pvalue <-
    function(x, digits = getOption("digits"), ...)
{
    frmt <- function(x, tol) {
        idx <- x < tol
        ndd <- n_decimal_digits(x)
        ## ndd is > 0 unless p-value = 0 (or 1)
        x[idx] <- if (ndd > 0) max(tol, 10^-ndd) else tol
        ## flag p-values below tolerance
        x <- format(x, digits = digits)
        x[idx] <- paste0("<", x[idx])
        x
    }
    ## set tolerance to either machine or resampling precision
    tol <- attr(x, "nresample") # non-NULL for approximate p-values
    tol <- if (is.null(tol)) eps else 1 / tol
    ## remove attributes
    mostattributes(x) <- list(
        dim = dim(x), dimnames = dimnames(x), names = names(x)
    )
    ## formatting using 'frmt()'
    if (!anyNA(x)) { # midpvalue() returns NA in some cases
        x[] <- if (length(dim(x)) > 0) # [] is used to *always* keep dimnames
                   ## handle each column separately just like 'print.default()'
                   apply(x, 2, frmt, tol = tol)
               else
                   frmt(x, tol = tol)
    }
    x
}

print.pvalue <-
    function (x, digits = getOption("digits"), quote = FALSE, right = TRUE, ...)
{
    ## print p-value...
    print(format(x, digits = digits), quote = quote, right = right, ...)
    ## ...and its confidence interval
    conf.int <- attr(x, "conf.int")
    if (!is.null(conf.int)) {
        conf.int <- list(conf.int = conf.int)
        class(conf.int) <- "ci"
        print(conf.int, digits = digits, ...)
    }
    invisible(x)
}

print.cutpoint <- function(x, ...) {
    cat(paste0("  ", dQuote("best"), " cutpoint: ", x$label, "\n"))
    if (hasName(x, "covariable"))
        cat(paste0("       covariable: ", x$covariable, "\n"))
}
