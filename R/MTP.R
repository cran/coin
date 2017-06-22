### single step maxT multiple testing procedure
singlestep <- function(object, ...) {

    ## reorder test statistics to ensure consistency with "global"/"step-down"
    switch(object@statistic@alternative,
           "two.sided" = {
               ts <- abs(statistic(object, type = "standardized"))
               o <- order(ts, decreasing = TRUE)}, # abs. largest ts first
           "greater" = {
               ts <- statistic(object, type = "standardized")
               o <- order(ts, decreasing = TRUE)}, # largest ts first
           "less" = {
               ts <- statistic(object, type = "standardized")
               o <- order(ts) # smallest ts first
           })

    ## iterate over unique test statistics only and remap
    pq <- length(ts)
    ots <- ts[o]
    idx <- c(which(ots[-1L] %NE% ots[-pq]), pq) # unique ts
    ret <- object@distribution@pvalue(ots[idx], ...)

    matrix(rep.int(ret, diff(c(0L, idx)))[order(o)], # remapping
           nrow = nrow(ts), ncol = ncol(ts), dimnames = dimnames(ts))
}

### algorithm 2.8 (Free Step-Down Resampling Method) in
### Westfall & Young (1993), page 66 _using standardized
### statistics instead of p-values_!
rsdmaxT <- function(pls, ts) {

    ## reorder simulations using (increasing) test statistics
    o <- order(ts) # smallest ts first
    pls <- pls[, o, drop = FALSE]

    ## algorithm 2.8 (Free Step-Down Resampling Method) in
    ## Westfall & Young (1993), page 66 _using standardized
    ## statistics instead of p-values_!
    if (ncol(pls) > 1) {
        for (j in 2:ncol(pls))
            pls[, j] <- pmax.int(pls[, j], pls[, j - 1])
    }
    ret <- rowMeans(t(pls) %GE% ts[o])
    for (i in (length(ret) - 1):1)
        ret[i] <- max(ret[i], ret[i + 1]) # enforce monotonicity, page 67

    matrix(ret[order(o)], nrow = nrow(ts), ncol = ncol(ts),
           dimnames = dimnames(ts))
}

### step-down using the asymptotic distribution
asdmaxT <- function(object) {

    ## reorder upper and/or lower limits using test statistics
    switch(object@statistic@alternative,
           "two.sided" = {
               ts <- abs(statistic(object, type = "standardized"))
               pq <- length(ts)
               o <- order(ts, decreasing = TRUE) # abs. largest ts first
               upper <- ts[o]
               lower <- -upper},
           "greater" = {
               ts <- statistic(object, type = "standardized")
               pq <- length(ts)
               o <- order(ts, decreasing = TRUE) # largest ts first
               upper <- ts[o]
               lower <- rep.int(-Inf, pq)},
           "less" = {
               ts <- statistic(object, type = "standardized")
               pq <- length(ts)
               o <- order(ts) # smallest ts first
               upper <- rep.int(Inf, pq)
               lower <- ts[o]})

    ## correlation matrix
    corr <- cov2cor(covariance(object))

    ## step-down based on multivariate normality
    ret <- numeric(pq)
    ret[1] <- pmvn(lower = lower[1], upper = upper[1],
                   mean = rep.int(0, pq), corr = corr,
                   conf.int = FALSE)
    if (pq > 1) {
        oo <- o
        for (i in 2:pq) {
            j <- rank(oo)[1] # reindexing needed in each step
            corr <- corr[-j, -j]
            oo <- oo[-1]
            ret[i] <- min(ret[i - 1],
                          pmvn(lower = lower[i], upper = upper[i],
                               mean = rep.int(0, length(oo)), corr = corr,
                               conf.int = FALSE))
        }
    }

    matrix(1 - ret[order(o)], nrow = nrow(ts), ncol = ncol(ts),
           dimnames = dimnames(ts))
}

### stepdown maxT multiple testing procedure
stepdown <- function(object, ...) {

    if (!(extends(class(object), "MaxTypeIndependenceTest")))
        stop(sQuote("object"), " is not of class ",
             sQuote("MaxTypeIndependenceTest"))

    if (extends(class(object@distribution), "AsymptNullDistribution"))
        asdmaxT(object)
    else {
        ## raw simulation results, scores have been handled already
        pls <- support(object, raw = TRUE)

        ## standardize
        dcov <- sqrt(variance(object))
        expect <- expectation(object)
        switch(object@statistic@alternative,
               "two.sided" = {
                   pls <- abs(t((pls - expect) / dcov))
                   ts <- abs(statistic(object, type = "standardized"))},
               "greater" = {
                   pls <- t((pls - expect) / dcov)
                   ts <- statistic(object, type = "standardized")},
               "less" = {
                   pls <- -t((pls - expect) / dcov)
                   ts <- -(statistic(object, type = "standardized"))})

        rsdmaxT(pls, ts)
    }
}

### Adjusted marginal p-values taking discreteness into account in the
### permutation case (Westfall and Wolfinger, 1997)
marginal <- function(object, bonferroni, stepdown, ...) {

    ## <FIXME> this should be possible when the _exact_ marginal
    ## distributions are available
    ## </FIXME>

    if (!(extends(class(object), "MaxTypeIndependenceTest")))
        stop(sQuote("object"), " is not of class ",
             sQuote("MaxTypeIndependenceTest"))

    if (extends(class(object@distribution), "AsymptNullDistribution")) {
        ## unadjusted p-values
        ts <- statistic(object, type = "standardized")
        ret <- switch(object@statistic@alternative,
                      "two.sided" = 2 * pmin.int(pnorm(ts), 1 - pnorm(ts)),
                      "greater"   = 1 - pnorm(ts),
                      "less"      = pnorm(ts))

        ## adjustment
        ret <- if (!stepdown) {
            if (bonferroni) pmin.int(1, length(ret) * ret) # Bonferroni
            else 1 - (1 - ret)^length(ret) # Sidak
        } else {
            n <- length(ret)
            o <- order(ret)
            if (bonferroni) # Bonferroni-Holm
                pmin.int(1, cummax((n - seq_len(n) + 1L) * ret[o])[order(o)])
            else # Sidak-Holm
                cummax(1 - (1 - ret[o])^(n - seq_len(n) + 1L))[order(o)]
        }

        matrix(ret, nrow = nrow(ts), ncol = ncol(ts), dimnames = dimnames(ts))
    } else {
        ## raw simulation results, scores have been handled already
        pls <- support(object, raw = TRUE)

        ## standardize
        dcov <- sqrt(variance(object))
        expect <- expectation(object)
        switch(object@statistic@alternative,
               "two.sided" = {
                   pls <- abs(t((pls - expect) / dcov))
                   ts <- abs(statistic(object, type = "standardized"))},
               "greater" = {
                   pls <- t((pls - expect) / dcov)
                   ts <- statistic(object, type = "standardized")},
               "less" = {
                   pls <- -t((pls - expect) / dcov)
                   ts <- -(statistic(object, type = "standardized"))})

        ## reorder simulations using the (decreasing) test statistics
        o <- order(ts, decreasing = TRUE) # largest ts first
        pls <- pls[, o, drop = FALSE]

        ## unadjusted p-values
        pu <- rowMeans(t(pls) %GE% ts[o])

        ## permutation distribution
        foo <- function(x, t) mean(x %GE% t)
        p <- vector(mode = "list", length = ncol(pls))
        for (i in 1:ncol(pls)) {
            ux <- unique(pls[, i])
            p[[i]] <- vapply(ux, foo, NA_real_, x = pls[, i])
        }

        ## discrete adjustment
        ret <- rep.int(1 - bonferroni, length(ts)) # zeros (ones) for Bonferroni (Sidak)
        for (i in 1:length(pu)) {
            qq <- if (stepdown) i else 1 # 'i' => successively smaller subsets
            for (q in qq:length(p)) {
                x <- p[[q]][p[[q]] <= pu[i]] # below eq. 2
                if (length(x) > 0) {
                    ret[i] <- if (bonferroni) ret[i] + max(x) # eq. 4
                              else ret[i] * (1 - max(x)) # eq. 2
                }
            }
        }
        ret <- if (!bonferroni) pmin.int(1 - ret, 1) else pmin.int(ret, 1)
        for (i in 2:length(ret))
            ret[i] <- max(ret[i - 1], ret[i]) # enforce monotonicity

        matrix(ret[order(o)], nrow = nrow(ts), ncol = ncol(ts),
               dimnames = dimnames(ts))
    }
}

### compute p-values under subset pivotality (Westfall, 1997)
npmcp <- function(object) {

    ## extract from object
    y <- object@statistic@y[[1]]
    x <- object@statistic@x[[1]]
    ytrafo <- object@statistic@ytrafo
    alternative <- object@statistic@alternative

    ## <FIXME> it is currently hard to ask a distribution object
    ## for its type (and arguments). Its a design bug.
    distribution <- object@call$distribution
    ## </FIXME>
    stand_tstat <- statistic(object, type = "standardized")
    tstat <- switch(alternative,
                    "less" = stand_tstat,
                    "greater" = -stand_tstat,
                    "two.sided" = -abs(stand_tstat))

    ## get contrast matrix from xtrans
    C <- attr(object@statistic@xtrans, "contrast")
    stopifnot(inherits(C, "matrix"))

    ## order test statistics, most "extreme" one comes first
    Corder <- C[order(tstat), , drop = FALSE]

    ## compute allowed subsets of hypotheses
    ## returns list consisting of lists (one for each rejection step of H0)
    ms <- multcomp:::maxsets(Corder)

    ## make sure 'object' isn't serialized along with 'foo'
    ## (otherwise parallel operation using snow clusters will be very slow)
    object <- NULL # was rm(object)
    ## alternatively we could pass all relevant objects to 'foo' and then
    ## associate it with the global environment instead:
    ## foo <- function(s, y, x, ytrafo, distribution, alternative) { ... }
    ## environment(foo) <- .GlobalEnv
    ## or simply define 'foo' out of 'npmcp'

    foo <- function(s) {
        Ctmp <- Corder[s, , drop = FALSE] # current allowed subset
        ## x levels in current subset
        xlev <- apply(Ctmp, MARGIN = 2, function(col) any(col != 0))

        it <- independence_test(y ~ x,
                                subset = x %in% names(xlev)[xlev], # relevant data subset
                                xtrafo = mcp_trafo(x = Ctmp),
                                ytrafo = ytrafo,
                                distribution = distribution,
                                alternative = alternative)
        pvalue(it)
    }

    p <- vapply(ms, function(sub) # for every list of allowed subsets
        max(vapply(sub, foo, NA_real_)), NA_real_) # for every subset

    for (i in 2:length(p))
        p[i] <- max(p[i-1], p[i]) # forces pvalue monotonicity

    matrix(p[rank(tstat)], dimnames = dimnames(tstat))
}

### unadjusted p-values
unadjusted <- function(object, ...) {

    if (extends(class(object@distribution), "AsymptNullDistribution")) {
        ts <- statistic(object, type = "standardized")
        ret <- switch(object@statistic@alternative,
                      "two.sided" = 2 * pmin.int(pnorm(ts), 1 - pnorm(ts)),
                      "greater"   = 1 - pnorm(ts),
                      "less"      = pnorm(ts))

        matrix(ret, nrow = nrow(ts), ncol = ncol(ts), dimnames = dimnames(ts))
    } else {
        ## raw simulation results, scores have been handled already
        pls <- support(object, raw = TRUE)

        ## standardize
        dcov <- sqrt(variance(object))
        expect <- expectation(object)
        switch(object@statistic@alternative,
               "two.sided" = {
                   pls <- abs((pls - expect) / dcov)
                   ts <- abs(statistic(object, type = "standardized"))},
               "greater" = {
                   pls <- (pls - expect) / dcov
                   ts <- statistic(object, type = "standardized")},
               "less" = {
                   pls <- -(pls - expect) / dcov
                   ts <- -(statistic(object, type = "standardized"))})

        ## unadjusted p-values
        matrix(rowMeans(pls %GE% as.vector(ts)),
               nrow = nrow(ts), ncol = ncol(ts), dimnames = dimnames(ts))
    }
}
