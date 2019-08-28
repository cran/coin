setGeneric(".confint",
    function(object1, object2, ...) {
        standardGeneric(".confint")
    }
)

setMethod(".confint",
    signature = list("ScalarIndependenceTestStatistic", "NullDistribution"),
    definition = function(object1, object2, parm, level, ...) {
        ## <FIXME> drop unused levels!
        if (!is_2sample(object1))
            warning(sQuote("object1"), " does not represent a two-sample problem")
        ## </FIXME>
        if (nlevels(object1@block) != 1L || !is_unity(object1@weights))
            stop("cannot compute confidence interval with blocks or weights")
        if (!(length(level) == 1L && level > 0 && level < 1))
            stop("level must be a single number between 0 and 1")

        parm <- match.arg(parm, choices = c("location", "scale"))
        location <- parm == "location"
        alpha <- 1 - level

        scores <- object1@y[[1L]]
        groups <- object1@xtrans[, 1L]
        ytrafo <- object1@ytrafo
        mu <- expectation(object1)
        sigma <- sqrt(variance(object1))

        ## raw data
        x <- sort(scores[groups > 0])
        y <- sort(scores[groups < 1])

        foo <- if (location) function(x, d) x - d
               else          function(x, d) x / d

        ## explicitly compute all possible steps
        steps <- outer(x, y, foo)
        if (!location)
            steps <- steps[steps >= 0]
        steps <- sort(unique(steps))
        steps <- steps[c(steps[-1L] %NE% steps[-length(steps)], TRUE)] # unique +/- eps

        ## computes the statistic under the alternative 'd'
        fs <- function(d)
            sum(ytrafo(data.frame(c(foo(x, d), y)))[seq_along(x)])

        ## we need to compute the statistics just to the right of
        ## each step
        jumps <- vapply(steps + min(diff(steps)) / 2, fs, NA_real_)

        ## determine if the statistics are in- or decreasing
        ## jumpsdiffs <- diff(jumps)
        increasing <- all(diff(jumps[c(1L, length(jumps))]) > 0)
        decreasing <- all(diff(jumps[c(1L, length(jumps))]) < 0)
        ## this is safe
        if (!(increasing || decreasing))
            stop("cannot compute confidence interval: ",
                 "the step function is not monotone")

        cci <- function(alpha) {
            ## the quantiles: reject iff
            ##   STATISTIC <  qlower OR
            ##   STATISTIC >= qupper
            qlower <- qperm(object2,     alpha / 2) * sigma + mu
            qupper <- qperm(object2, 1 - alpha / 2) * sigma + mu
            ## Check if the statistic exceeds both quantiles first.
            if (qlower < min(jumps) || qupper > max(jumps)) {
                warning("cannot compute confidence interval")
                return(c(NA, NA))
            }

            if (increasing) {
                ## do NOT reject for all steps with
                ##   STATISTICS >= qlower AND
                ##   STATISTICS <  qupper
                ## but the open right interval ends with the
                ## step with STATISTIC == qupper
                c(min(steps[jumps %GE% qlower]), min(steps[jumps %GT% qupper]))
            } else {
                ## do NOT reject for all steps with
                ##   STATISTICS >= qlower AND
                ##   STATISTICS <  qupper
                ## but the open left interval ends with the
                ## step with STATISTIC == qupper
                c(min(steps[jumps %LE% qupper]), min(steps[jumps %LT% qlower]))
            }
        }

        cint <- switch(object1@alternative,
                    "two.sided" = cci(alpha),
                    "greater"   = c(cci(alpha * 2)[1L], Inf),
                    "less"      = c(if (location) -Inf else 0, cci(alpha * 2)[2L])
                )
        attr(cint, "conf.level") <- level

        ## was: median(steps) which will not work for blocks etc.
        ESTIMATE <- if (increasing)
                        c(min(steps[jumps %GE% mu]), min(steps[jumps %GT% mu]))
                    else
                        c(min(steps[jumps %LE% mu]), min(steps[jumps %LT% mu]))
        ESTIMATE <- mean(ESTIMATE, na.rm = TRUE)
        names(ESTIMATE) <- if (location) "difference in location"
                           else          "ratio of scales"

        list(conf.int = cint, estimate = ESTIMATE)
    }
)

setMethod(".confint",
    signature = list("ScalarIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, parm, level, ...) {
        ## <FIXME> drop unused levels!
        if (!is_2sample(object1))
            warning(sQuote("object1"), " does not represent a two-sample problem")
        ## </FIXME>
        if (nlevels(object1@block) != 1L || !is_unity(object1@weights))
            stop("cannot compute confidence interval with blocks or weights")
        if (!(length(level) == 1L && level > 0 && level < 1))
            stop("level must be a single number between 0 and 1")

        parm <- match.arg(parm, choices = c("location", "scale"))
        location <- parm == "location"
        alpha <- 1 - level

        scores <- object1@y[[1L]]
        groups <- object1@xtrans[, 1L]
        ytrafo <- object1@ytrafo
        mu <- expectation(object1)
        sigma <- sqrt(variance(object1))

        ## raw data
        x <- sort(scores[groups > 0])
        y <- sort(scores[groups < 1])

        foo <- if (location) function(x, d) x - d
               else          function(x, d) x / d

        ## approximate the steps
        ## Here we search the root of the function 'fs' on the set
        ## c(mumin, mumax).
        ##
        ## This returns a value from c(mumin, mumax) for which
        ## the standardized statistic is equal to the
        ## quantile zq.  This means that the statistic is not
        ## within the critical region, and that implies that '
        ## is a confidence limit for the median.

        fs <- function(d, zq)
            (sum(ytrafo(data.frame(c(foo(x, d), y)))[seq_along(x)]) - mu) / sigma - zq

        if (location) {
            mumin <- min(x) - max(y)
            mumax <- max(x) - min(y)
        } else {
            srangepos <- NULL
            srangeneg <- NULL
            if (any(xGT0 <- x %GT% 0) && any(yGT0 <- y %GT% 0)) {
                xGT0 <- x[xGT0]
                yGT0 <- y[yGT0]
                srangepos <-
                    c(min(xGT0, na.rm = TRUE) / max(yGT0, na.rm = TRUE),
                      max(xGT0, na.rm = TRUE) / min(yGT0, na.rm = TRUE))
            }
            if (any(xLE0 <- x %LE% 0) && any(yLT0 <- y %LT% 0)) {
                xLE0 <- x[xLE0]
                yLT0 <- y[yLT0]
                srangeneg <-
                    c(min(xLE0, na.rm = TRUE) / max(yLT0, na.rm = TRUE),
                      max(xLE0, na.rm = TRUE) / min(yLT0, na.rm = TRUE))
            }
            if (any(is.infinite(c(srangepos, srangeneg))))
                stop("cannot compute asymptotic confidence set or estimator")
            srange <- c(srangepos, srangeneg)
            mumin <- min(srange)
            mumax <- max(srange)
        }

        cci <- function(alpha) {
            ## Check if the statistic exceeds both quantiles
            ## first: otherwise 'uniroot' won't work anyway
            statu <- fs(mumin, zq = qperm(object2,     alpha / 2))
            statl <- fs(mumax, zq = qperm(object2, 1 - alpha / 2))
            if (sign(statu) == sign(statl)) {
                warning("samples differ in location: ",
                        "cannot compute confidence set, returning NA")
                return(c(NA, NA))
            }
            u <- uniroot(fs, c(mumin, mumax),
                         zq = qperm(object2,     alpha / 2), tol = sqrt_eps)$root
            l <- uniroot(fs, c(mumin, mumax),
                         zq = qperm(object2, 1 - alpha / 2), tol = sqrt_eps)$root
            ## The process of the statistics does not need to be
            ## increasing: sort is ok here.
            sort(c(u, l))
        }

        cint <- switch(object1@alternative,
                    "two.sided" = cci(alpha),
                    "greater"   = c(cci(alpha * 2)[1L], Inf),
                    "less"      = c(if (location) -Inf else 0, cci(alpha * 2)[2L])
                )
        attr(cint, "conf.level") <- level

        ## Check if the statistic exceeds both quantiles first.
        statu <- fs(mumin, zq = 0)
        statl <- fs(mumax, zq = 0)
        ESTIMATE <- if (sign(statu) == sign(statl)) {
                        warning("cannot compute estimate, returning NA")
                        NA
                    } else
                        uniroot(fs, c(mumin, mumax), zq = 0, tol = sqrt_eps)$root
        names(ESTIMATE) <- if (location) "difference in location"
                           else          "ratio of scales"

        list(conf.int = cint, estimate = ESTIMATE)
    }
)

setMethod(".confint",
    signature = list("ScalarIndependenceTest", "missing"),
    definition = function(object1, object2, parm, level, ...) {
        callGeneric(object1@statistic, object1@distribution, parm, level, ...)
    }
)


### CI for a binomial parameter
confint_binom <-
    function(x, n, level = 0.95, method = c("exact", "mid-p"), tol = eps)
{
    method <- match.arg(method)

    RET <- if (x >= 0 && x <= n) {
               alpha <- 1 - level
               if (method == "exact") {
                   ## exact Clopper-Pearson interval
                   c(if (x == 0) 0 else qbeta(    alpha / 2, x    , n - x + 1),
                     if (x == n) 1 else qbeta(1 - alpha / 2, x + 1, n - x    ))
               } else {
                   ## mid-p interval (see Berry and Armitage, 1995)
                   if (x == 0) {
                       c(0, 1 - alpha^(1 / n))
                   } else if (x == n) {
                       c(alpha^(1 / n), 1)
                   } else {
                       f <- function(p, a)
                           ## 0.5 * dbinom(...) + pbinom(..., lower.tail = FALSE)
                           mean(pbinom(c(x, x - 1), n, p, lower.tail = FALSE)) - a
                       c(uniroot(f, c(0, 1), a =     alpha / 2, tol = tol)$root,
                         uniroot(f, c(0, 1), a = 1 - alpha / 2, tol = tol)$root)
                   }
               }
           } else {
               stop(sQuote("x"), " must be larger or equal to 0 and",
                    " smaller or equal to ", sQuote("n"))
           }
    attr(RET, "conf.level") <- level
    RET
}


###
### Currently unused
###
simconfint_location <- function(object, level = 0.95,
    approx = FALSE, ...) {

    if (!(is_Ksample(object@statistic) &&
        inherits(object, "MaxTypeIndependenceTest")))
        stop(sQuote("object"), " is not of class ",
             dQuote("MaxTypeIndependenceTest"),
             " representing a K-sample problem")

    xtrans <- object@statistic@xtrans
    if (!all(apply(xtrans, 2L, function(x) all(x %in% c(-1, 0, 1)))))
        stop("only differences are allowed as contrasts")

    estimate <- c()
    lower <- c()
    upper <- c()

    ## transform max(abs(x))-type distribution into a
    ## distribution symmetric around zero
    nnd <- object@distribution
    nnd@q <- function(p) {
        pp <- p
        if (p > 0.5)
            pp <- 1 - pp
        q <- qperm(object@distribution, 1 - pp)
        if (p < 0.5)
            -q
        else
            q
    }

    for (i in seq_len(ncol(xtrans))) {
        thisset <- abs(xtrans[, i]) > 0
        ip <- new("IndependenceProblem",
                  object@statistic@x[thisset, , drop = FALSE],
                  object@statistic@y[thisset, , drop = FALSE],
                  object@statistic@block[thisset])

        it <- independence_test(ip, teststat = "scalar",
            distribution = "none", alternative = "two.sided",
            ytrafo = object@statistic@ytrafo, ...)

        ci <- .confint(it@statistic, nnd, parm = "location", level = level, ...)
        estimate <- c(estimate, ci$estimate)
        lower <- c(lower, ci$conf.int[1L])
        upper <- c(upper, ci$conf.int[2L])
    }
    RET <- data.frame(Estimate = estimate, lower = lower, upper = upper)
    colnames(RET)[2L:3L] <-
        paste(c((1 - level) / 2, 1 - (1 - level) / 2) * 100, "%")
    rownames(RET) <- colnames(object@statistic@xtrans)
    attr(RET, "conf.level") <- level
    RET
}
