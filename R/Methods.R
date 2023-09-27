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
        p <- function(q) pnorm(q)                        # implicitly vectorized
        q <- function(p) qnorm(p)                        # implicitly vectorized
        d <- function(x) dnorm(x)                        # implicitly vectorized
        pvalue <- function(q) {                          # implicitly vectorized
            switch(object@alternative,
                "less"      = p(q),
                "greater"   = 1 - p(q),
                "two.sided" = 2 * pmin(p(q), 1 - p(q))
            )
        }

        new("AsymptNullDistribution",
            seed = NA_integer_,
            p = p,
            q = q,
            d = d,
            pvalue = pvalue,
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            size = function(alpha, type) NA,
            support = function() NA,
            name = "Univariate Normal Distribution")
    }
)

### method for max-type test statistics
setMethod("AsymptNullDistribution",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(object, ...) {
        seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (is.null(seed)) {
            runif(1L)
            seed <- .GlobalEnv[[".Random.seed"]]
        }

        corr <- cov2cor(covariance(object, partial = FALSE))
        pq <- nrow(corr)

        p_fun <- function(q, conf.int, ...) {
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
        q_fun <- function(p, ...) {
            qmvn(p, mean = rep.int(0, pq), corr = corr, ...)
        }

        p <- function(q, ...) {
            setAttributes(vapply(q, p_fun, NA_real_, conf.int = FALSE, ...,
                                 USE.NAMES = FALSE),
                          attributes(q))
        }
        q <- function(p, ...) {
            setAttributes(vapply(p, q_fun, NA_real_, ..., USE.NAMES = FALSE),
                          attributes(p))
        }
        pvalue <- function(q, ...) {
            if (length(q) < 2L) {
                RET <- 1 - p_fun(q, conf.int = TRUE, ...)
                ci <- 1 - attr(RET, "conf.int")[2L:1L]
                attr(ci, "conf.level") <-
                    attr(attr(RET, "conf.int"), "conf.level")
                attr(RET, "conf.int") <- ci
                RET
            } else
                1 - vapply(q, p_fun, NA_real_, conf.int = FALSE, ...)
        }

        new("AsymptNullDistribution",
            seed = seed,
            p = p,
            q = q,
            d = function(x) NA,
            pvalue = pvalue,
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            size = function(alpha, type) NA,
            support = function() NA,
            name = "Multivariate Normal Distribution",
            parameters = list(corr = corr))
    }
)

### method for quad-type test statistics
setMethod("AsymptNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, ...) {
        df <- object@df
        p <- function(q) pchisq(q, df = df)              # implicitly vectorized
        q <- function(p) qchisq(p, df = df)              # implicitly vectorized
        d <- function(x) dchisq(x, df = df)              # implicitly vectorized
        pvalue <- function(q) 1 - p(q)                   # implicitly vectorized

        new("AsymptNullDistribution",
            seed = NA_integer_,
            p = p,
            q = q,
            d = d,
            pvalue = pvalue,
            midpvalue = function(q) NA,
            pvalueinterval = function(q) NA,
            size = function(alpha, type) NA,
            support = function() NA,
            name = "Chi-Squared Distribution",
            parameters = list(df = df))
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
    definition = function(object, nresample = 10000L, B, ...) {
        ## <DEPRECATED>
        if (!missing(B)) {
            warning(sQuote("B"), " is deprecated; use ", sQuote("nresample"),
                    " instead")
            nresample <- B
        }
        ## </DEPRECATED>
        seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (is.null(seed)) {
            runif(1L)
            seed <- .GlobalEnv[[".Random.seed"]]
        }

        pls <- MonteCarlo(object@xtrans, object@ytrans, object@block,
                          object@weights, nresample, ...)
        pls <- (pls - as.vector(.expectation(object, partial = FALSE))) /
                   sqrt(as.vector(.variance(object, partial = FALSE)))

        p_fun <- function(q) {
            mean(pls %LE% q)
        }
        d_fun <- function(x) {
            mean(pls %EQ% x)
        }
        pvalue_fun <- function(q, conf.int) {
            RET <- switch(object@alternative,
                       "less"      = mean(pls %LE% q),
                       "greater"   = mean(pls %GE% q),
                       "two.sided" = mean(abs(pls) %GE% abs(q))
                   )
            if (conf.int) {
                attr(RET, "conf.int") <-
                    confint_binom(round(RET * nresample), nresample,
                                  level = 0.99, method = "exact")
            }
            RET
        }
        midpvalue_fun <- function(q, conf.int, z) {
            RET <- pvalue_fun(q, conf.int = FALSE) - z *
                     if (object@alternative == "two.sided")
                         d_fun(-q) + d_fun(q) # both tails
                     else
                         d_fun(q)
            if (conf.int) {
                attr(RET, "conf.int") <-
                    confint_binom(round(RET * nresample), nresample,
                                  level = 0.99, method = "mid-p")
            }
            RET
        }

        p <- function(q) {
            setAttributes(vapply(q, p_fun, NA_real_, USE.NAMES = FALSE),
                          attributes(q))
        }
        q <- function(p) {                               # implicitly vectorized
            setAttributes(quantile(pls, probs = p, names = FALSE, type = 1L),
                          attributes(p))
        }
        d <- function(x) {
            setAttributes(vapply(x, d_fun, NA_real_, USE.NAMES = FALSE),
                          attributes(x))
        }
        pvalue <- function(q) {
            if (length(q) < 2L)
                pvalue_fun(q, conf.int = TRUE)
            else
                vapply(q, pvalue_fun, NA_real_, conf.int = FALSE)
        }
        midpvalue <- function(q) {
            if (length(q) < 2L)
                midpvalue_fun(q, conf.int = TRUE, z = 0.5)
            else
                vapply(q, midpvalue_fun, NA_real_, conf.int = FALSE, z = 0.5)
        }
        pvalueinterval <- function(q) {
            if (length(q) < 2L)
                midpvalue_fun(q, conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
            else
                vapply(q, midpvalue_fun, c(NA_real_, NA_real_),
                       conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
        }
        support <- function(raw = FALSE) {
            if (raw)
                pls
            else {
                ## NOTE: '%NE%' is expensive, so drop duplicates first
                pls <- sort(unique(pls))
                pls[c(pls[-1L] %NE% pls[-length(pls)], TRUE)] # unique +/- eps
            }
        }
        size <- function(alpha, type) {
            spt <- support()
            pv <- if (type == "mid-p-value")
                      midpvalue(spt)
                  else
                      pvalue(spt)
            vapply(alpha, function(a) sum(d(spt[pv %LE% a])), NA_real_)
        }

        new("ApproxNullDistribution",
            seed = seed,
            nresample = nresample,
            p = p,
            q = q,
            d = d,
            pvalue = pvalue,
            midpvalue = midpvalue,
            pvalueinterval = pvalueinterval,
            size = size,
            support = support,
            name = "Monte Carlo Distribution")
    }
)

### method for max-type test statistics
setMethod("ApproxNullDistribution",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(object, nresample = 10000L, B, ...) {
        ## <DEPRECATED>
        if (!missing(B)) {
            warning(sQuote("B"), " is deprecated; use ", sQuote("nresample"),
                    " instead")
            nresample <- B
        }
        ## </DEPRECATED>
        seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (is.null(seed)) {
            runif(1L)
            seed <- .GlobalEnv[[".Random.seed"]]
        }

        pls <- MonteCarlo(object@xtrans, object@ytrans, object@block,
                          object@weights, nresample, ...)
        pls <- (pls - as.vector(.expectation(object, partial = FALSE))) /
                   sqrt(as.vector(.variance(object, partial = FALSE)))

        mpls <- switch(object@alternative,
                    "less"      = colMins(pls),
                    "greater"   = colMaxs(pls),
                    "two.sided" = colMaxs(abs(pls))
                )

        p_fun <- function(q) {
            switch(object@alternative,
                "less"      = mean(mpls %GE% q),
                "greater"   = mean(mpls %LE% q),
                "two.sided" = mean(mpls %LE% q)
            )
        }
        d_fun <- function(x) {
            mean(mpls %EQ% x)
        }
        pvalue_fun <- function(q, conf.int) {
            RET <- switch(object@alternative,
                       "less"      = mean(mpls %LE% q),
                       "greater"   = mean(mpls %GE% q),
                       "two.sided" = mean(mpls %GE% q)
                   )
            if (conf.int) {
                attr(RET, "conf.int") <-
                    confint_binom(round(RET * nresample), nresample,
                                  level = 0.99, method = "exact")
            }
            RET
        }
        midpvalue_fun <- function(q, conf.int, z) {
            RET <- pvalue_fun(q, conf.int = FALSE) - z *
                     if (object@alternative == "two.sided")
                         d_fun(-q) + d_fun(q) # both tails
                     else
                         d_fun(q)
            if (conf.int) {
                attr(RET, "conf.int") <-
                    confint_binom(round(RET * nresample), nresample,
                                  level = 0.99, method = "mid-p")
            }
            RET
        }

        p <- function(q) {
            setAttributes(vapply(q, p_fun, NA_real_, USE.NAMES = FALSE),
                          attributes(q))
        }
        q <- function(p) {                               # implicitly vectorized
            setAttributes(quantile(mpls, probs = p, names = FALSE, type = 1L),
                          attributes(p))
        }
        d <- function(x) {
            setAttributes(vapply(x, d_fun, NA_real_, USE.NAMES = FALSE),
                          attributes(x))
        }
        pvalue <- function(q) {
            if (length(q) < 2L)
                pvalue_fun(q, conf.int = TRUE)
            else
                vapply(q, pvalue_fun, NA_real_, conf.int = FALSE)
        }
        midpvalue <- function(q) {
            if (length(q) < 2L)
                midpvalue_fun(q, conf.int = TRUE, z = 0.5)
            else
                vapply(q, midpvalue_fun, NA_real_, conf.int = FALSE, z = 0.5)
        }
        pvalueinterval <- function(q) {
            if (length(q) < 2L)
                midpvalue_fun(q, conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
            else
                vapply(q, midpvalue_fun, c(NA_real_, NA_real_),
                       conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
        }
        support <- function(raw = FALSE) {
            if (raw)
                pls
            else {
                ## NOTE: '%NE%' is expensive, so drop duplicates first
                mpls <- sort(unique(mpls))
                mpls[c(mpls[-1L] %NE% mpls[-length(mpls)], TRUE)] # unique +/- eps
            }
        }
        size <- function(alpha, type) {
            spt <- support()
            pv <- if (type == "mid-p-value")
                      midpvalue(spt)
                  else
                      pvalue(spt)
            vapply(alpha, function(a) sum(d(spt[pv %LE% a])), NA_real_)
        }

        new("ApproxNullDistribution",
            seed = seed,
            nresample = nresample,
            p = p,
            q = q,
            d = d,
            pvalue = pvalue,
            midpvalue = midpvalue,
            pvalueinterval = pvalueinterval,
            size = size,
            support = support,
            name = "Monte Carlo Distribution")
    }
)

### method for quad-type test statistics
setMethod("ApproxNullDistribution",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, nresample = 10000L, B, ...) {
        ## <DEPRECATED>
        if (!missing(B)) {
            warning(sQuote("B"), " is deprecated; use ", sQuote("nresample"),
                    " instead")
            nresample <- B
        }
        ## </DEPRECATED>
        seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (is.null(seed)) {
            runif(1L)
            seed <- .GlobalEnv[[".Random.seed"]]
        }

        pls <- MonteCarlo(object@xtrans, object@ytrans, object@block,
                          object@weights, nresample, ...)
        pls <- .Call(R_quadform,
                     linstat = pls,
                     expect = .expectation(object, partial = FALSE),
                     MPinv_sym = object@covarianceplus)

        p_fun <- function(q) {
            mean(pls %LE% q)
        }
        d_fun <- function(x) {
            mean(pls %EQ% x)
        }
        pvalue_fun <- function(q, conf.int) {
            RET <- mean(pls %GE% q)
            if (conf.int) {
                attr(RET, "conf.int") <-
                    confint_binom(round(RET * nresample), nresample,
                                  level = 0.99, method = "exact")
            }
            RET
        }
        midpvalue_fun <- function(q, conf.int, z) {
            RET <- pvalue_fun(q, conf.int = FALSE) - z * d_fun(q)
            if (conf.int) {
                attr(RET, "conf.int") <-
                    confint_binom(round(RET * nresample), nresample,
                                  level = 0.99, method = "mid-p")
            }
            RET
        }

        p <- function(q) {
            setAttributes(vapply(q, p_fun, NA_real_, USE.NAMES = FALSE),
                          attributes(q))
        }
        q <- function(p) {                               # implicitly vectorized
            setAttributes(quantile(pls, probs = p, names = FALSE, type = 1L),
                          attributes(p))
        }
        d <- function(x) {
            setAttributes(vapply(x, d_fun, NA_real_, USE.NAMES = FALSE),
                          attributes(x))
        }
        pvalue <- function(q) {
            if (length(q) < 2L)
                pvalue_fun(q, conf.int = TRUE)
            else
                vapply(q, pvalue_fun, NA_real_, conf.int = FALSE)
        }
        midpvalue <- function(q) {
            if (length(q) < 2L)
                midpvalue_fun(q, conf.int = TRUE, z = 0.5)
            else
                vapply(q, midpvalue_fun, NA_real_, conf.int = FALSE, z = 0.5)
        }
        pvalueinterval <- function(q) {
            if (length(q) < 2L)
                midpvalue_fun(q, conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
            else
                vapply(q, midpvalue_fun, c(NA_real_, NA_real_),
                       conf.int = FALSE, z = c("p_0" = 1, "p_1" = 0))
        }
        support <- function(raw = FALSE) {
            if (raw)
                pls
            else {
                ## NOTE: '%NE%' is expensive, so drop duplicates first
                pls <- sort(unique(pls))
                pls[c(pls[-1L] %NE% pls[-length(pls)], TRUE)] # unique +/- eps
            }
        }
        size <- function(alpha, type) {
            spt <- support()
            pv <- if (type == "mid-p-value")
                      midpvalue(spt)
                  else
                      pvalue(spt)
            vapply(alpha, function(a) sum(d(spt[pv %LE% a])), NA_real_)
        }

        new("ApproxNullDistribution",
            seed = seed,
            nresample = nresample,
            p = p,
            q = q,
            d = d,
            pvalue = pvalue,
            midpvalue = midpvalue,
            pvalueinterval = pvalueinterval,
            size = size,
            support = support,
            name = "Monte Carlo Distribution")
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
            if (NCOL(object@ytrans) > 1L)
                stop("cannot compute exact distribution with multivariate scores")
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


### method for extraction of confidence intervals
setMethod("confint",
    signature = "IndependenceTest",
    definition = function(object, parm, level = 0.95, ...) {
        stop("cannot compute confidence interval for objects of class ",
             dQuote(class(object)))
    }
)

setMethod("confint",
    signature = "ScalarIndependenceTestConfint",
    definition = function(object, parm, level = 0.95, ...) {
        ci <- if (missing(level))
                  object@confint(object@conf.level)
              else
                  object@confint(level)
        class(ci) <- "ci"
        ci
    }
)
