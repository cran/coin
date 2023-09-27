### Pearson's chi-squared test
chisq_test <- function(object, ...) UseMethod("chisq_test")

chisq_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("chisq_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

chisq_test.table <- function(object, ...) {

    do.call(chisq_test,
            c(object = table2IndependenceProblem(object), list(...)))
}

chisq_test.IndependenceProblem <- function(object, ...) {

    check <- function(object) {
        if (!is_contingency(object))
            stop(sQuote("object"),
                 " does not represent a contingency problem")
        if (nlevels(object@block) != 1)
            stop(sQuote("object"), " contains blocks: use ",
                 sQuote("cmh_test"), " instead")
        TRUE
    }
    n <- sum(object@weights)

    args <- setup_args()
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quadratic"
    ## distribution must be checked
    args$distribution <- check_distribution_arg(args$distribution)
    ## alternative is needed later
    args$alternative <- match.arg(args$alternative,
                                  c("two.sided", "less", "greater"))

    ## transform data if requested and setup a test problem
    object <- new("IndependenceTestProblem", object, args$xtrafo, args$ytrafo)

    if (!check(object))
        stop(sQuote("check"), " failed")

    ## use the classical chisq statistic based on Pearson
    ## residuals (O - E)^2 / E
    ## see Th. 3.1 and its proof in Strasser & Weber (1999).
    object <- new("IndependenceLinearStatistic", object)
    object@covariance <- object@covariance * (n - 1) / n
    object <-
        if (args$teststat == "scalar") {
            object <-
                new("ScalarIndependenceTestStatistic", object, args$alternative)
            new("ScalarIndependenceTest", statistic = object,
                distribution = args$distribution(object))
        } else {
            if (args$alternative != "two.sided")
                warning(sQuote("alternative"),
                        " is ignored for quadratic test statistics")
            object <- new("QuadTypeIndependenceTestStatistic", object)
            new("QuadTypeIndependenceTest", statistic = object,
                distribution = args$distribution(object))
        }

    if (is_doubly_ordered(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else if (is_singly_ordered(object@statistic))
        object@method <- "Generalized Pearson Chi-Squared Test"
    else
        object@method <- "Pearson Chi-Squared Test"

    object@call <- match.call()

    object
}


### generalized Cochran-Mantel-Haenzel test
cmh_test <- function(object, ...) UseMethod("cmh_test")

cmh_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("cmh_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

cmh_test.table <- function(object, ...) {

    do.call(cmh_test,
            c(object = table2IndependenceProblem(object), list(...)))
}

cmh_test.IndependenceProblem <- function(object, ...) {

    args <- setup_args(
        check = function(object) {
            if (!is_contingency(object))
                stop(sQuote("object"),
                     " does not represent a contingency problem")
            TRUE
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quadratic"

    object <- do.call(independence_test, c(object = object, args))

    if (is_doubly_ordered(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else
        object@method <- "Generalized Cochran-Mantel-Haenszel Test"

    object
}


### linear-by-linear association test
lbl_test <- function(object, ...) UseMethod("lbl_test")

lbl_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("lbl_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

lbl_test.table <- function(object, ...) {

    do.call(lbl_test,
            c(object = table2IndependenceProblem(object), list(...)))
}

lbl_test.IndependenceProblem <- function(object, ...) {

    ## convert factors to ordered
    object@x[] <- lapply(object@x, function(x)
        if (is.factor(x) && nlevels(x) > 2) as.ordered(x) else x)
    object@y[] <- lapply(object@y, function(y)
        if (is.factor(y) && nlevels(y) > 2) as.ordered(y) else y)

    args <- setup_args(
        teststat = "scalar",
        check = function(object) {
            if (!is_doubly_ordered(object))
                stop(sQuote("object"),
                     " does not represent a problem with ordered data")
            TRUE
        }
    )

    object <- do.call(independence_test, c(object = object, args))

    object@method <- "Linear-by-Linear Association Test"

    object
}
