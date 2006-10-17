
ft <- function(test, formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    d <- formula2data(formula, data, subset, weights = weights, ...)
    ip <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl,
              weights = d$w)
    RET <- do.call(test, c(list(object = ip), list(...)))
    return(RET)
}

### a generic test procedure for classical (and not so classical) tests
independence_test <- function(object, ...) UseMethod("independence_test")

independence_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("independence_test", formula, data, subset, weights, ...)
}

independence_test.table <- function(object, 
    distribution = c("asymptotic", "approximate"), ...) {

    distribution <- check_distribution_arg(distribution, 
                                           c("asymptotic", "approximate"))
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (class(distribution) == "asymptotic") {
        df <- as.data.frame(object)
        if (ncol(df) == 3)
            ip <- new("IndependenceProblem", x = df[1], y = df[2], block = NULL, 
                      weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]], weights = df[["Freq"]])
        }
    } else {
        df <- table2df(object)
        if (ncol(df) == 3) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]])
        } else {
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = NULL) 
        }
    }
    ### </FIXME>
    RET <- do.call("independence_test", 
                   c(list(object = ip, distribution = distribution), 
                   list(...)))
    return(RET)
}

independence_test.IndependenceProblem <- function(object,
    teststat = c("max", "quad", "scalar"),
    distribution = c("asymptotic", "approximate", "exact"), 
    alternative = c("two.sided", "less", "greater"), 
    xtrafo = trafo, ytrafo = trafo, scores = NULL, check = NULL, ...) {

    addargs <- list(...)
    if (length(addargs) > 0) 
        warning("additional arguments ", 
                paste(names(addargs), collapse = ", "),
                " will be ignored")

    ### just for backward compatibility
    teststat <- match.arg(teststat, choices = c("maxtype", "quadtype", "scalar"), 
                          several.ok = TRUE)
    if (teststat[1] == "maxtype") teststat <- "max"
    if (teststat[1] == "quadtype") teststat <- "quad"
    teststat <- match.arg(teststat)
    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    ### convert factors to ordered and attach scores if requested
    object <- setscores(object, scores)

    ### transform data if requested and setup a test problem
    itp <- new("IndependenceTestProblem", object, xtrafo = xtrafo, 
               ytrafo = ytrafo, ...)

    if (!is.null(check)) {
        if (is.function(check)) {
            if (!check(itp))
                stop(sQuote("check"), " failed")
        } else {
            stop(sQuote("check"), " is not a function")
        }
    }

    ### check type of test statistic and alternative
    scalar <- is_scalar(itp)

    if (!scalar) {
        if (teststat == "scalar") {
            warning("Length linear statistic > 1, using ",
                    sQuote("max"), "-type test statistic")
            teststat <- "max"
        }
    } else {
        if (teststat == "max") teststat <- "scalar"
    }
    if (alternative != "two.sided" && teststat == "quad")
        warning(sQuote("alternative"), " is ignored for ", 
                teststat, " type test statistics")

    ### compute linear statistic, conditional expectation and
    ### conditional covariance
    its <- new("IndependenceTestStatistic", itp, 
        varonly = class(distribution) == "approximate" && teststat == "max")

    ### compute test statistic and corresponding null distribution
    RET <- switch(teststat,
        "scalar" = {
            ts <- new("ScalarIndependenceTestStatistic", its, 
                      alternative = alternative)

            nd <- switch(class(distribution),
                "asymptotic" = do.call("AsymptNullDistribution", 
                                       c(list(object = ts), distribution)),
                "exact"  = do.call("ExactNullDistribution", 
                                   c(list(object = ts), distribution)),
                "approximate" = do.call("ApproxNullDistribution", 
                                        c(list(object = ts), distribution))
            )
            new("ScalarIndependenceTest", statistic = ts, distribution = nd)
        },
        "max" = {
            ts <- new("MaxTypeIndependenceTestStatistic", its, 
                      alternative = alternative)
            nd <- switch(class(distribution),
                "asymptotic" = do.call("AsymptNullDistribution", 
                                       c(list(object = ts), distribution)),
                "exact"  = do.call("ExactNullDistribution", 
                                   c(list(object = ts), distribution)),
                "approximate" = do.call("ApproxNullDistribution", 
                                        c(list(object = ts), distribution))
                      )
            new("MaxTypeIndependenceTest", statistic = ts, distribution = nd)
        },
        "quad" = {
            ts <- new("QuadTypeIndependenceTestStatistic", its)
            nd <- switch(class(distribution),
                "asymptotic" = do.call("AsymptNullDistribution", 
                                       c(list(object = ts), distribution)),
                "exact"  = do.call("ExactNullDistribution", 
                                   c(list(object = ts), distribution)),
                "approximate" = do.call("ApproxNullDistribution", 
                                        c(list(object = ts), distribution))
            )
            new("QuadTypeIndependenceTest", statistic = ts, 
                distribution = nd)
        })

    RET@method <- "General Independence Test"

    ### return object inheriting from class `IndependenceTest'
    return(RET)
}


### OK, OK, here is the most prominent one ...
wilcox_test <- function(object, ...) UseMethod("wilcox_test")

wilcox_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("wilcox_test", formula, data, subset, weights, ...)
}

wilcox_test.IndependenceProblem <- function(object,  
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"), 
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_2sample(object) && is_numeric_y(object)))
            stop(sQuote("object"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    RET <- independence_test(object, teststat = "scalar", 
        alternative = alternative, distribution = distribution, 
        ytrafo = function(data) trafo(data, numeric_trafo = rank), 
        check = check, ...)

    RET@nullvalue <- 0
    RET@method <- "Wilcoxon Mann-Whitney Rank Sum Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_location(RET@statistic, RET@distribution, 
                             level = level, 
                             approx = (class(distribution) == "asymptotic"))
        RET@conf.level <- conf.level
    }

   return(RET)
}


### normal quantiles (van der Waerden) test
normal_test <- function(object, ...) UseMethod("normal_test")

normal_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("normal_test", formula, data, subset, weights, ...)
}   

normal_test.IndependenceProblem <- function(object,  
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"), 
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_2sample(object) && is_numeric_y(object)))
            stop(sQuote("object"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    RET <- independence_test(object, teststat = "scalar", 
        alternative = alternative, distribution = distribution, 
        ytrafo = function(data) trafo(data, numeric_trafo = function(x)
            normal_trafo(x, ties.method = ties.method)), 
        check = check, ...)

    RET@nullvalue <- 0
    RET@method <- "Normal Quantile (van der Waerden) Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_location(RET@statistic, RET@distribution,
                             level = level, 
                             approx = (class(distribution) == "asympt"))
        RET@conf.level <- conf.level
    }
    return(RET)
}


### median test
median_test <- function(object, ...) UseMethod("median_test")

median_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("median_test", formula, data, subset, weights, ...)
}   

median_test.IndependenceProblem <- function(object,     
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"), 
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_2sample(object) && is_numeric_y(object)))
            stop(sQuote("object"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    RET <- independence_test(object, teststat = "scalar", 
        alternative = alternative, distribution = distribution, 
        ytrafo = function(data) trafo(data, numeric_trafo = median_trafo), 
        check = check, ...)
 
    RET@nullvalue <- 0
    RET@method <- "Median Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_location(RET@statistic, RET@distribution,
                             level = level, 
                             approx = (class(distribution) == "asymptotic"))
        RET@conf.level <- conf.level
    }
    return(RET)
}


### Ansari-Bradley test
ansari_test <- function(object, ...) UseMethod("ansari_test")

ansari_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("ansari_test", formula, data, subset, weights, ...)
}   

ansari_test.IndependenceProblem <- function(object,
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"), 
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {     

    check <- function(object) {
        if (!(is_2sample(object) && is_numeric_y(object)))
            stop(sQuote("object"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    if (alternative == "less") {
        alternative <- "greater"
    } else {
        if (alternative == "greater") 
            alternative <- "less"
    }
    distribution <- check_distribution_arg(distribution)

    RET <- independence_test(object, teststat = "scalar",
        alternative = alternative, distribution = distribution,
        ytrafo = function(data) trafo(data, numeric_trafo = function(x)
            ansari_trafo(x, ties.method = ties.method)), 
        check = check, ...)
 
    RET@nullvalue <- 1
    RET@method <- "Ansari-Bradley Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_scale(RET@statistic, RET@distribution,
                          level = level, 
                          approx = (class(distribution) == "asymptotic"))
        RET@conf.level <- conf.level
    }
    return(RET)
}


### Survival tests -> Logrank only, for the moment
surv_test <- function(object, ...) UseMethod("surv_test")

surv_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("surv_test", formula, data, subset, weights, ...)
}
    
surv_test.IndependenceProblem <- function(object,  
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"), 
    ties.method = c("logrank", "HL"), ...) {

    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)
    ties.method <- match.arg(ties.method)

    ytrafo <- function(data) trafo(data, surv_trafo = function(x)
        logrank_trafo(x, ties.method = ties.method))

    check <- function(object) {
        if (!(is_Ksample(object) && is_censored_y(object)))
            stop(sQuote("object"), 
                 " does not represent a K sample problem with censored data")
        return(TRUE)
    }

    scalar <- FALSE
    if (is.factor(object@x[[1]])) scalar <- nlevels(object@x[[1]]) == 2

    RET <- independence_test(object, 
        teststat = ifelse(scalar, "scalar", "quad"), 
        distribution = distribution, check = check, ytrafo = ytrafo, 
        ...)
 
    if (extends(class(RET@statistic), "ScalarIndependenceTest"))
        RET@nullvalue <- 0

    if (is_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association (Tarone-Ware) Test"
    else
        RET@method <- "Logrank Test"
    return(RET)
}


### Kruskal-Wallis test
kruskal_test <- function(object, ...) UseMethod("kruskal_test")

kruskal_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("kruskal_test", formula, data, subset, weights, ...)
}   

kruskal_test.IndependenceProblem <- function(object,  
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"), " does not represent a K sample problem")
        return(TRUE)
    }
 
    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    RET <- independence_test(object, 
        distribution = distribution, teststat = "quad",
        ytrafo = function(data) trafo(data, numeric_trafo = rank), 
        check = check, ...)

    if (is_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else
        RET@method <- "Kruskal-Wallis Test"
    return(RET)
}


### Fligner test
fligner_test <- function(object, ...) UseMethod("fligner_test")

fligner_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("fligner_test", formula, data, subset, weights, ...)
}   

fligner_test.IndependenceProblem <- function(object,  
    ties.method = c("mid-ranks", "average-scores"),
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"), " does not represent a K sample problem")
        if (is_ordered(object))
            stop(colnames(object@x), " is an ordered factor")
        return(TRUE)
    }
 
    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    ### eliminate location differences (see `stats/R/fligner.test')
    object@y[[1]] <- object@y[[1]] - 
        tapply(object@y[[1]], object@x[[1]], median)[object@x[[1]]]

    RET <- independence_test(object,  
        distribution = distribution, teststat = "quad",
        ytrafo = function(data) trafo(data, numeric_trafo = function(x)
            fligner_trafo(x, ties.method = ties.method)), 
        check = check, ...)

    RET@method <- "Fligner-Killeen Test"
    return(RET)
}


### Spearman test
spearman_test <- function(object, ...) UseMethod("spearman_test")

spearman_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("spearman_test", formula, data, subset, weights, ...)
}   

spearman_test.IndependenceProblem <- function(object, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_corr(object))
            stop(sQuote("object"), 
                 " does not represent a univariate correlation problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    RET <- independence_test(object, 
        teststat = "scalar", alternative = alternative, 
        distribution = distribution, 
        xtrafo = function(data) trafo(data, numeric_trafo = rank),
        ytrafo = function(data) trafo(data, numeric_trafo = rank), 
        check = check, ...)

    RET@nullvalue <- 0
    RET@method <- "Spearman Correlation Test"
    return(RET)
}


### Generalised Cochran-Mantel-Haenzel Test
cmh_test <- function(object, ...) UseMethod("cmh_test")

cmh_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("cmh_test", formula, data, subset, weights, ...)
}   

cmh_test.table <- function(object, distribution = c("asymptotic", "approximate"), ...) {

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (class(distribution) == "asymptotic") {
        df <- as.data.frame(object)
        if (ncol(df) == 3)
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = NULL, weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]], weights = df[["Freq"]])
        }
    } else {
        df <- table2df(object)
        if (ncol(df) == 3) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]])
        } else {
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = NULL) 
        }
    }
    ### </FIXME>
    RET <- do.call("cmh_test", c(list(object = ip, distribution = distribution), 
                   list(...)))
    return(RET)
}

cmh_test.IndependenceProblem <- function(object, 
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_contingency(object))
            stop(sQuote("object"), " does not represent a contingency problem")
        return(TRUE)
    }

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    RET <- independence_test(object, 
        teststat = "quad", distribution = distribution, check = check, 
        ...)

    if (is_ordered(RET@statistic)) 
        RET@method <- "Linear-by-Linear Association Test"
    else
        RET@method <- "Generalised Cochran-Mantel-Haenszel Test"
    return(RET)
}


### Pearsons Chi-Squared Test
chisq_test <- function(object, ...) UseMethod("chisq_test")

chisq_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("chisq_test", formula, data, subset, weights, ...)
}   

chisq_test.table <- function(object, distribution = c("asymptotic", "approximate"), ...) {

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (class(distribution) == "asymptotic") {
        df <- as.data.frame(object)
        if (ncol(df) == 3)
            ip <- new("IndependenceProblem", x = df[1], y = df[2],
                      block = NULL, weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2],
                      block = df[[3]], weights = df[["Freq"]])
        }
    } else {
        df <- table2df(object)
        if (ncol(df) == 3) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2],
                      block = df[[3]])
        } else {
            ip <- new("IndependenceProblem", x = df[1], y = df[2],
                      block = NULL)
        }
    }
    ### </FIXME>
    RET <- do.call("chisq_test", c(list(object = ip, distribution = distribution), 
                   list(...)))
    return(RET)
}

chisq_test.IndependenceProblem <- function(object,  
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_contingency(object))
            stop(sQuote("object"), " does not represent a contingency problem")
        if (nlevels(object@block) != 1)
            stop(sQuote("object"), " contains blocks: use ", 
                 sQuote("cmh_test"), " instead")
        return(TRUE)
    }
    n <- sum(object@weights)

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    RET <- independence_test(object, 
        teststat = "quad", distribution = "asymptotic", check = check, ...)

    ### use the classical chisq statistic based on Pearson 
    ### residuals (O - E)^2 / E
    ### see Th. 3.1 and its proof in Strasser & Weber (1999).

    RET@statistic@teststatistic <- 
        RET@statistic@teststatistic * n / (n - 1)
    RET@statistic@covariance <- 
        new("CovarianceMatrix", covariance(RET) * (n - 1) / n)
    RET@statistic@covarianceplus <- MPinv(covariance(RET))$MPinv

    if (class(distribution) == "approximate") {
        nd <- do.call("ApproxNullDistribution", 
                      c(list(object = RET@statistic), distribution))
        RET <- new("QuadTypeIndependenceTest", statistic = RET@statistic,
                distribution = nd)
    }

    if (is_ordered(RET@statistic)) 
        RET@method <- "Linear-by-Linear Association Test"
    else
        RET@method <- "Pearson's Chi-Squared Test"
    return(RET)
}


### Linear-by-Linear Association Test
lbl_test <- function(object, ...) UseMethod("lbl_test")

lbl_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("lbl_test", formula, data, subset, weights, ...)
}   

lbl_test.table <- function(object, 
    distribution = c("asymptotic", "approximate"), ...) {

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (class(distribution) == "asymptotic") {
        df <- as.data.frame(object)
        if (nlevels(df[[1]]) > 2) df[[1]] <- ordered(df[[1]])
        if (nlevels(df[[2]]) > 2) df[[2]] <- ordered(df[[2]])
        if (ncol(df) == 3)
            ip <- new("IndependenceProblem", x = df[1], y = df[2], block = NULL, 
                      weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]], weights = df[["Freq"]])
        }
    } else {
        df <- table2df(object)
        if (nlevels(df[[1]]) > 2) df[[1]] <- ordered(df[[1]])
        if (nlevels(df[[2]]) > 2) df[[2]] <- ordered(df[[2]])
        if (ncol(df) == 3) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]])
        } else {
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = NULL) 
        }
    }
    RET <- do.call("lbl_test", c(list(object = ip, distribution = distribution), 
                   list(...)))
    return(RET)
}


lbl_test.IndependenceProblem <- function(object, 
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_ordered(object))
            stop(sQuote("object"), 
                 " does not represent a problem with ordered data")
        return(TRUE)
    }

    ### convert factors to ordered
    object@x <- as.data.frame(lapply(object@x, 
        function(x) if (is.factor(x) && nlevels(x) > 2) {
                        return(ordered(x))
                    } else {
                        return(x)
                    }))
    object@y <- as.data.frame(lapply(object@y, 
        function(x) if (is.factor(x) && nlevels(x) > 2) {
                        return(ordered(x))
                    } else {
                        return(x)
                    }))

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    RET <- do.call("independence_test", 
        c(list(object = object, teststat = "quad", distribution = distribution, 
               check = check), list(...)))

    RET@method <- "Linear-by-Linear Association Test"
    return(RET)
}


### permutation test without transformations
oneway_test <- function(object, ...) UseMethod("oneway_test")

oneway_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("oneway_test", formula, data, subset, weights, ...)
}

oneway_test.IndependenceProblem <- function(object, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asymptotic", "approximate", "exact"), ...) {

    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"), " does not represent a K sample problem")
        return(TRUE)
    }

    RET <- independence_test(object, 
        alternative = alternative, distribution = distribution, 
        check = check, ...)

    if (is_scalar(RET@statistic))
        RET@nullvalue <- 0
    RET@method <- paste(ifelse(length(table(object@x[[1]])) == 2, "2-", "K-"), 
                        paste("Sample Permutation Test"), sep = "")
    return(RET)
}


### Contrast test
contrast_test <- function(object, ...) UseMethod("contrast_test")

contrast_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("contrast_test", formula, data, subset, weights, ...)
}

contrast_test.IndependenceProblem <- function(object, 
    cmatrix, distribution = c("asymptotic", "approximate"), ...) {

    if (!(ncol(object@x) == 1 && is.factor(object@x[[1]])))
        stop(sQuote("object@x"), " is not univariate or a factor")

    if  (!is.matrix(cmatrix) || nrow(cmatrix) != nlevels(object@x[[1]]))
        stop(sQuote("cmatrix"), " is not a matrix with ", 
             nlevels(object@x), " rows")

    if (is.null(colnames(cmatrix)))
        colnames(cmatrix) <- paste("C", 1:ncol(cmatrix), sep = "")

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    xtrafo <- function(data) trafo(data) %*% cmatrix

    RET <- independence_test(object, teststat = "max",
        distribution = distribution, xtrafo = xtrafo, ...)
    RET@method <- "General Contrast Test"
    
    return(RET)
}


### Maxstat test
maxstat_test <- function(object, ...) UseMethod("maxstat_test")

maxstat_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("maxstat_test", formula, data, subset, weights, ...)
}

maxstat_test.IndependenceProblem <- function(object, 
    distribution = c("asymptotic", "approximate"), 
    teststat = c("max", "quad"), 
    minprob = 0.1, maxprob = 0.9, ...) {

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))
    teststat <- match.arg(teststat)

    ORDERED <- sapply(object@x, is.ordered)
    lev <- lapply(object@x, levels)
    for (i in which(ORDERED)) class(object@x[[i]]) <- "numeric"

    mm <- function(x) maxstat_trafo(x, minprob = minprob, maxprob = maxprob)
    fmm <- function(x) fmaxstat_trafo(x, minprob = minprob, maxprob = maxprob)
    xtrafo <- function(data) trafo(data, numeric_trafo = mm, factor_trafo = fmm)

    RET <- independence_test(object, teststat = teststat,
        distribution = distribution, xtrafo = xtrafo, ...)

    ### estimate cutpoint
    wm <- which.max(apply(abs(statistic(RET, "standardized")), 1, max))
    whichvar <- attr(RET@statistic@xtrans, "assign")[wm]
    maxcontr <- RET@statistic@xtrans[,wm]
    if (is.factor(object@x[[whichvar]])) {
        estimate <- levels(object@x[[whichvar]][maxcontr > 0][, drop = TRUE])
    } else {
        estimate <- max(object@x[[whichvar]][maxcontr > 0])
        if (ORDERED[whichvar]) estimate <- lev[[whichvar]][estimate]
    }
    if (ncol(object@x) > 1) {
        estimate <- list(covariable = colnames(object@x)[whichvar], cutpoint = estimate)
    } else {
        estimate <- list(cutpoint = estimate)
    }
    RET@statistic@estimates <- list(estimate = estimate)
    RET@method <- "Maxstat Test"
    
    return(RET)
}


### a generic test procedure for classical (and not so classical) tests
symmetry_test <- function(object, ...) UseMethod("symmetry_test")

symmetry_test.formula <- function(formula, data = list(), subset = NULL,
    ...) {

    d <- formula2data(formula, data, subset, ...)
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("symmetry_test", c(list(object = sp), list(...)))
    return(RET)

}

symmetry_test.SymmetryProblem <- function(object,
    teststat = c("max", "quad", "scalar"),
    distribution = c("asymptotic", "approximate", "exact"), 
    alternative = c("two.sided", "less", "greater"), 
    xtrafo = trafo, ytrafo = trafo, scores = NULL, 
    check = NULL, ...) {
    class(object) <- "IndependenceProblem"
    independence_test(object, teststat, distribution, alternative, xtrafo,
                      ytrafo, scores, check, ...)
}

symmetry_test.table <- function(object, ...) {
    df <- table2df_sym(object)
    sp <- new("SymmetryProblem", x = df["groups"], y = df["response"])
    RET <- do.call("symmetry_test", c(list(object = sp), list(...))) 
    return(RET)
}

### Friedman-Test
friedman_test <- function(object, ...) UseMethod("friedman_test")

friedman_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("friedman_test", c(list(object = sp), list(...)))
    return(RET)
}   

friedman_test.SymmetryProblem <- function(object, 
    distribution = c("asymptotic", "approximate"), ...) {
    
    if (!is_completeblock(object))
        stop("Not an unreplicated complete block design")

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    RET <- symmetry_test(object, 
        distribution = distribution, teststat = "quad", 
        ytrafo = function(data) 
            trafo(data, numeric_trafo = rank, block = object@block), 
        ...)

    if (is_ordered(RET@statistic))  
        RET@method <- "Page Test"
    else
        RET@method <- "Friedman Test"
    return(RET)
}

### Marginal-Homogeneity-Test
mh_test <- function(object, ...) UseMethod("mh_test")

mh_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("mh_test", c(list(object = sp), list(...)))
    return(RET)
}   

mh_test.table <- function(object, ...) {
    df <- table2df_sym(object)
    sp <- new("SymmetryProblem", x = df["groups"], y = df["response"])
    RET <- do.call("mh_test", c(list(object = sp), list(...))) 
    return(RET)
}

mh_test.SymmetryProblem <- function(object, 
    distribution = c("asymptotic", "approximate"), ...) {
    
    if (!is_completeblock(object))
        stop("Not an unreplicated complete block design")
    if (ncol(object@y) != 1 || !is.factor(object@y[[1]]))
        stop("Response variable is not a factor")

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate"))

    addargs <- list(...)
    scores <- addargs$scores
    if (!is.null(addargs$scores)) {
        if (length(scores) > 1)
            stop("length of ", sQuote("scores"), " must be equal one")
        names(scores) <- "response"
        addargs$scores <- NULL
    }
 
    RET <- do.call("symmetry_test", 
        c(list(object = object, distribution = distribution, 
               teststat = "quad", scores = scores), addargs))

    if (is_ordered(RET@statistic))  
        RET@method <- "Marginal-Homogeneity Test for Ordered Data"
    else
        RET@method <- "Marginal-Homogeneity Test"
    return(RET)
}


### Wilcoxon-Signed-Rank Test
wilcoxsign_test <- function(object, ...) UseMethod("wilcoxsign_test")

wilcoxsign_test.formula <- function(formula, data = list(), 
                                    subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    ip <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("wilcoxsign_test", c(list(object = ip), list(...)))
    return(RET)
}   

wilcoxsign_test.IndependenceProblem <- function(object, 
    distribution = c("asymptotic", "approximate", "exact"), ...) {

    y <- object@y[[1]]
    x <- object@x[[1]]
    if (!is.numeric(x))
        stop(sQuote("x"), " is not a numeric variable")
    if (!is.numeric(y))
        stop(sQuote("y"), " is not a numeric variable")

    block <- gl(length(x), 2)
    diffs <- x - y
    diffs <- diffs[abs(diffs) > 0]
    block <- gl(length(diffs), 2)
    pos <- rank(abs(diffs)) * (diffs > 0)
    neg <- rank(abs(diffs)) * (diffs < 0)
    yy <- drop(as.vector(t(cbind(pos, neg))))
    xx <- factor(rep(c("pos", "neg"), length(diffs)))

    distribution <- check_distribution_arg(distribution, 
        values = c("asymptotic", "approximate", "exact"))

    ip <- new("IndependenceProblem", x = data.frame(x = xx), 
              y = data.frame(y = yy), block = block)

    RET <- independence_test(ip,
        distribution = distribution, teststat = "scalar", ...)

    RET@method <- "Wilcoxon-Signed-Rank Test"
    RET@nullvalue <- 0
    return(RET)
}
