### generic method for asymptotic null distributions
setGeneric("AsymptNullDistribution",
    function(object, ...) {
        standardGeneric("AsymptNullDistribution")
    }
)

### method for scalar test statistics
setMethod("AsymptNullDistribution",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(object, ...) {
        p <- function(q) pnorm(q)
        q <- function(p) qnorm(p)
        d <- function(x) dnorm(x)
        pvalue <- function(q) {
            switch(object@alternative,
                "less"      = p(q),
                "greater"   = 1 - p(q),
                "two.sided" = 2 * pmin.int(p(q), 1 - p(q))
            )
        }

        new("AsymptNullDistribution",
            seed = NA_integer_,
            p = p,                           # implicitly vectorized
            q = q,                           # implicitly vectorized
            d = d,                           # implicitly vectorized
            pvalue = pvalue,                 # implicitly vectorized
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            support = function() NA,
            name = "Univariate Normal Distribution")
    }
)

### method for max-type test statistics
setMethod("AsymptNullDistribution",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(object, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        corr <- cov2cor(covariance(object))
        pq <- length(expectation(object))

        p <- function(q, conf.int, ...) {
            switch(object@alternative,
                "less"      = pmvn(lower = q, upper = Inf,
                                   mean = rep.int(0, pq),
                                   corr = corr, conf.int = conf.int, ...),
                "greater"   = pmvn(lower = -Inf, upper = q,
                                   mean = rep.int(0, pq),
                                   corr = corr, conf.int = conf.int, ...),
                "two.sided" = pmvn(lower = -abs(q), upper = abs(q),
                                   mean = rep.int(0, pq),
                                   corr = corr, conf.int = conf.int, ...)
            )
        }
        q <- function(p, ...) {
            qmvn(p, mean = rep.int(0, pq), corr = corr, ...)
        }
        pvalue <- function(q, conf.int, ...) {
            RET <- 1 - p(q, conf.int, ...)
            if (conf.int) {
                ci <- 1 - attr(RET, "conf.int")[2L:1L]
                attr(ci, "conf.level") <-
                    attr(attr(RET, "conf.int"), "conf.level")
                attr(RET, "conf.int") <- ci
                class(RET) <- "MCp"
            }
            RET
        }

        new("AsymptNullDistribution",
            seed = seed,
            p = function(q, ...) {
                vapply(q, p, NA_real_, conf.int = FALSE, ...)
            },
            q = function(p, ...) {
                vapply(p, q, NA_real_, ...)
            },
            d = function(x) NA,
            pvalue = function(q, ...) {
                if (length(q) < 2L)
                    pvalue(q, conf.int = TRUE, ...)
                else
                    vapply(q, pvalue, NA_real_, conf.int = FALSE, ...)
            },
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            support = function() NA,
            name = "Multivariate Normal Distribution",
            parameters = list(corr = corr))
    }
)

### method for quad-type test statistics
setMethod("AsymptNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, ...) {
        p <- function(q) pchisq(q, df = object@df)
        q <- function(p) qchisq(p, df = object@df)
        d <- function(d) dchisq(d, df = object@df)
        pvalue <- function(q) 1 - p(q)

        new("AsymptNullDistribution",
            seed = NA_integer_,
            p = p,                           # implicitly vectorized
            q = q,                           # implicitly vectorized
            d = d,                           # implicitly vectorized
            pvalue = pvalue,                 # implicitly vectorized
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            support = function() NA,
            name = "Chi-Squared Distribution",
            parameters = list(df = object@df))
    }
)


### generic method for exact null distributions
setGeneric("ExactNullDistribution",
    function(object, ...) {
        standardGeneric("ExactNullDistribution")
    }
)

### method for scalar test statistics
setMethod("ExactNullDistribution",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(object,
        algorithm = c("auto", "shift", "split-up"), ...) {
            algorithm <- match.arg(algorithm)
            if (object@paired) {
                if (algorithm == "split-up")
                    stop("split-up algorithm not implemented for paired samples")
                int <- is_integer(object@ytrans[, 1L], ...)
                if (int)
                    SR_shift_1sample(object, fact = attr(int, "fact"))
                else
                    stop("cannot compute exact distribution with real-valued scores")
            } else if (is_2sample(object)) {
                if (algorithm == "split-up")
                    vdW_split_up_2sample(object)
                else {
                    int <- is_integer(object@ytrans[, 1L], ...)
                    if (int)
                        SR_shift_2sample(object, fact = attr(int, "fact"))
                    else if (algorithm == "auto")
                        vdW_split_up_2sample(object)
                    else
                        stop("cannot compute exact distribution with real-valued scores")
                }
            } else
                stop(sQuote("object"), " is not a two-sample problem")
        }
)

### method for quad-type test statistics
setMethod("ExactNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object,
        algorithm = c("auto", "shift", "split-up"), ...) {
            algorithm <- match.arg(algorithm)
            if (object@paired) {
                if (algorithm == "split-up")
                    stop("split-up algorithm not implemented for paired samples")
                int <- is_integer(object@ytrans[, 1L], ...)
                if (int)
                    SR_shift_1sample(object, fact = attr(int, "fact"))
                else
                    stop("cannot compute exact distribution with real-valued scores")
            } else if (is_2sample(object)) {
                if (algorithm == "split-up")
                    stop("split-up algorithm not implemented for quadratic tests")
                else {
                    int <- is_integer(object@ytrans[, 1L], ...)
                    if (int)
                        SR_shift_2sample(object, fact = attr(int, "fact"))
                    else if (algorithm == "auto")
                        stop("split-up algorithm not implemented for quadratic tests")
                    else
                        stop("cannot compute exact distribution with real-valued scores")
                }
            } else
                stop(sQuote("object"), " is not a two-sample problem")
        }
)


### generic method for approximate null distributions
setGeneric("ApproxNullDistribution",
    function(object, ...) {
        standardGeneric("ApproxNullDistribution")
    }
)

### method for scalar test statistics
setMethod("ApproxNullDistribution",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(object, B = 10000, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        plsraw <-
            MonteCarlo(object@xtrans, object@ytrans, as.integer(object@block),
                       object@weights, as.integer(B), ...)

        ## <FIXME> can transform p, q, x instead of those </FIXME>
        pls <- sort((plsraw - expectation(object)) / sqrt(variance(object)))

        p <- function(q) {
            mean(pls %LE% q)
        }
        q <- function(p) {
            quantile(pls, probs = p, names = FALSE, type = 1L)
        }
        d <- function(x) {
            mean(pls %EQ% x)
        }
        pvalue <- function(q, conf.int) {
            RET <- switch(object@alternative,
                       "less"      = mean(pls %LE% q),
                       "greater"   = mean(pls %GE% q),
                       "two.sided" = mean(abs(pls) %GE% abs(q))
                   )
            if (conf.int) {
                attr(RET, "conf.int") <- confint_binom(round(RET * B), B)
                class(RET) <- "MCp"
            }
            RET
        }
        midpvalue <- function(q, conf.int, z) {
            RET <- pvalue(q, conf.int = FALSE) - z *
                     if (object@alternative == "two.sided")
                         d(-q) + d(q) # both tails
                     else
                         d(q)
            if (conf.int) {
                attr(RET, "conf.int") <- confint_midp(round(RET * B), B)
                class(RET) <- "MCp"
            }
            RET
        }

        new("ApproxNullDistribution",
            seed = seed,
            p = function(q) {
                vapply(q, p, NA_real_)
            },
            q = q,                           # implicitly vectorized
            d = function(x) {
                vapply(x, d, NA_real_)
            },
            pvalue = function(q) {
                if (length(q) < 2L)
                    pvalue(q, conf.int = TRUE)
                else
                    vapply(q, pvalue, NA_real_, conf.int = FALSE)
            },
            midpvalue = function(q) {
                if (length(q) < 2L)
                    midpvalue(q, conf.int = TRUE, z = 0.5)
                else
                    vapply(q, midpvalue, NA_real_, conf.int = FALSE, z = 0.5)
            },
            pvalueinterval = function(q) {
                if (length(q) < 2L)
                    midpvalue(q, conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
                else
                    vapply(q, midpvalue, c(NA_real_, NA_real_), conf.int = FALSE,
                           z = c("p_0" = 1, "p_1" = 0))
            },
            support = function(raw = FALSE) {
                if (raw)
                    plsraw
                else
                    pls[c(pls[-1L] %NE% pls[-length(pls)], TRUE)] # keep unique
            },
            name = "Monte Carlo Distribution")
    }
)

### method for max-type test statistics
setMethod("ApproxNullDistribution",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(object, B = 10000, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        plsraw <-
            MonteCarlo(object@xtrans, object@ytrans, as.integer(object@block),
                       object@weights, as.integer(B), ...)

        pls <- (plsraw - expectation(object)) / sqrt(variance(object))

        ## <FIXME>
        ## pls is a rather large object (potentially)
        ## try not to copy it too often -- abs() kills you
        ## </FIXME>

        pmaxmin <- function() {
            pls <- switch(object@alternative,
                       "less"      = do.call("pmin.int", as.data.frame(t(pls))),
                       "greater"   = do.call("pmax.int", as.data.frame(t(pls))),
                       "two.sided" = do.call("pmax.int", as.data.frame(t(abs(pls))))
                   )
            sort(pls)
        }

        p <- function(q) {
            switch(object@alternative,
                "less"      = mean(colSums(pls %GE% q) == nrow(pls)),
                "greater"   = mean(colSums(pls %LE% q) == nrow(pls)),
                "two.sided" = mean(colSums(abs(pls) %LE% q) == nrow(pls))
            )
        }
        q <- function(p) {
            quantile(pmaxmin(), probs = p, names = FALSE, type = 1L)
        }
        d <- function(x) {
            mean(pmaxmin() %EQ% x)
        }
        pvalue <- function(q, conf.int) {
            RET <- switch(object@alternative,
                       "less"      = mean(colSums(pls %LE% q) > 0),
                       "greater"   = mean(colSums(pls %GE% q) > 0),
                       "two.sided" = mean(colSums(abs(pls) %GE% q) > 0)
                   )
            if (conf.int) {
                attr(RET, "conf.int") <- confint_binom(round(RET * B), B)
                class(RET) <- "MCp"
            }
            RET
        }
        midpvalue <- function(q, conf.int, z) {
            RET <- pvalue(q, conf.int = FALSE) - z *
                     if (object@alternative == "two.sided")
                         d(-q) + d(q) # both tails
                     else
                         d(q)
            if (conf.int) {
                attr(RET, "conf.int") <- confint_midp(round(RET * B), B)
                class(RET) <- "MCp"
            }
            RET
        }

        new("ApproxNullDistribution",
            seed = seed,
            p = function(q) {
                vapply(q, p, NA_real_)
            },
            q = q,                           # implicitly vectorized
            d = function(x) {
                vapply(x, d, NA_real_)
            },
            pvalue = function(q) {
                if (length(q) < 2L)
                    pvalue(q, conf.int = TRUE)
                else
                    vapply(q, pvalue, NA_real_, conf.int = FALSE)
            },
            midpvalue = function(q) {
                if (length(q) < 2L)
                    midpvalue(q, conf.int = TRUE, z = 0.5)
                else
                    vapply(q, midpvalue, NA_real_, conf.int = FALSE, z = 0.5)
            },
            pvalueinterval = function(q) {
                if (length(q) < 2L)
                    midpvalue(q, conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
                else
                    vapply(q, midpvalue, c(NA_real_, NA_real_), conf.int = FALSE,
                           z = c("p_0" = 1, "p_1" = 0))
            },
            support = function(raw = FALSE) {
                if (raw)
                    plsraw
                else {
                    tmp <- pmaxmin()
                    tmp[c(tmp[-1L] %NE% tmp[-length(tmp)], TRUE)] # keep unique
                }
            },
            name = "Monte Carlo Distribution")
    }
)

### method for quad-type test statistics
setMethod("ApproxNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, B = 10000, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        plsraw <-
            MonteCarlo(object@xtrans, object@ytrans, as.integer(object@block),
                       object@weights, as.integer(B), ...)

        pls <- plsraw - expectation(object)
        pls <- sort(rowSums(crossprod(pls, object@covarianceplus) * t(pls)))

        p <- function(q) {
            mean(pls %LE% q)
        }
        q <- function(p) {
            quantile(pls, probs = p, names = FALSE, type = 1L)
        }
        d <- function(x) {
            mean(pls %EQ% x)
        }
        pvalue <- function(q, conf.int) {
            RET <- mean(pls %GE% q)
            if (conf.int) {
                attr(RET, "conf.int") <- confint_binom(round(RET * B), B)
                class(RET) <- "MCp"
            }
            RET
        }
        midpvalue <- function(q, conf.int, z) {
            RET <- pvalue(q, conf.int = FALSE) - z * d(q)
            if (conf.int) {
                attr(RET, "conf.int") <- confint_midp(round(RET * B), B)
                class(RET) <- "MCp"
            }
            RET
        }

        new("ApproxNullDistribution",
            seed = seed,
            p = function(q) {
                vapply(q, p, NA_real_)
            },
            q = q,                           # implicitly vectorized
            d = function(x) {
                vapply(x, d, NA_real_)
            },
            pvalue = function(q) {
                if (length(q) < 2L)
                    pvalue(q, conf.int = TRUE)
                else
                    vapply(q, pvalue, NA_real_, conf.int = FALSE)
            },
            midpvalue = function(q) {
                if (length(q) < 2L)
                    midpvalue(q, conf.int = TRUE, z = 0.5)
                else
                    vapply(q, midpvalue, NA_real_, conf.int = FALSE, z = 0.5)
            },
            pvalueinterval = function(q) {
                if (length(q) < 2L)
                    midpvalue(q, conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
                else
                    vapply(q, midpvalue, c(NA_real_, NA_real_), conf.int = FALSE,
                           z = c("p_0" = 1, "p_1" = 0))
            },
            support = function(raw = FALSE) {
                if (raw)
                    plsraw
                else
                    pls[c(pls[-1L] %NE% pls[-length(pls)], TRUE)] # keep unique
            },
            name = "Monte Carlo Distribution")
    }
)


### S3 method for extraction of confidence intervals
confint.ScalarIndependenceTestConfint <-
    function(object, parm, level = 0.95, ...) {
        x <- if ("level" %in% names(match.call()))
                 object@confint(level)
             else
                 object@confint(object@conf.level)
        class(x) <- "ci"
        x
}
