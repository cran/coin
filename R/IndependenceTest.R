### a generic test procedure for classical (and not so classical) tests
independence_test <- function(object, ...) UseMethod("independence_test")

independence_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("independence_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

independence_test.table <- function(object, ...) {

    do.call(independence_test,
            c(object = table2IndependenceProblem(object), list(...)))
}

independence_test.IndependenceProblem <- function(object,
    teststat = c("maximum", "quadratic", "scalar"),
    distribution = c("asymptotic", "approximate", "exact", "none"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL, check = NULL,
    ...) {

    if (...length() > 0L)
        warning("additional arguments ",
                paste0(sQuote(...names()), collapse = ", "),
                " will be ignored")

    teststat <- match.arg(teststat)
    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    ## convert factors to ordered and attach scores if requested
    if (!is.null(scores))
        object <- setscores(object, scores)

    ## transform data if requested and setup a test problem
    object <- new("IndependenceTestProblem", object, xtrafo = xtrafo,
                  ytrafo = ytrafo, ...)

    if (!is.null(check)) {
        if (is.function(check)) {
            if (!check(object))
                stop(sQuote("check"), " failed")
        } else {
            stop(sQuote("check"), " is not a function")
        }
    }

    ## check type of test statistic and alternative
    if (!is_scalar(object)) {
        if (teststat == "scalar") {
            warning("linear statistic has length > 1: ",
                    "using maximum test statistic instead")
            teststat <- "maximum"
        }
    } else {
        if (teststat == "maximum") teststat <- "scalar"
    }
    if (alternative != "two.sided" && teststat == "quadratic")
        warning(sQuote("alternative"),
                " is ignored for quadratic test statistic")

    ## compute linear statistic, conditional expectation and
    ## conditional covariance
    object <- new("IndependenceLinearStatistic", object)

    ## compute test statistic and corresponding null distribution
    ## return object inheriting from class "IndependenceTest"
    switch(teststat,
        "scalar" = {
            object <- new("ScalarIndependenceTestStatistic", object,
                          alternative = alternative, paired = FALSE)
            new("ScalarIndependenceTest", statistic = object,
                distribution = distribution(object), call = match.call())
        },
        "maximum" = {
            object <- new("MaxTypeIndependenceTestStatistic", object,
                          alternative = alternative)
            new("MaxTypeIndependenceTest", statistic = object,
                distribution = distribution(object), call = match.call())
        },
        "quadratic" = {
            object <- new("QuadTypeIndependenceTestStatistic", object,
                          paired = FALSE)
            new("QuadTypeIndependenceTest", statistic = object,
                distribution = distribution(object), call = match.call())
        }
    )
}
