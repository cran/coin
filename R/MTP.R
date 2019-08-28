### joint distribution-based max-T multiple testing procedures
### (Westfall and Young, 1993)
setGeneric("joint",
    function(object1, object2, ...) {
        standardGeneric("joint")
    }
)

setMethod("joint",
    signature = list("MaxTypeIndependenceTestStatistic", "NullDistribution"),
    definition = function(object1, object2, stepdown, ...) {
        if (!stepdown) {
            switch(object1@alternative,
                "less" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z)                    # smallest z first
                },
                "greater" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z, decreasing = TRUE) # largest z first
                },
                "two.sided" = {
                    z <- abs(statistic(object1, type = "standardized"))
                    o <- order(z, decreasing = TRUE) # abs. largest z first
                }
            )
            ## compute p-values for unique test statistics only and remap
            RET <- z[o]
            pq <- length(RET)
            idx <- c(which(RET[-1L] %NE% RET[-pq]), pq) # unique +/- eps
            RET <- rep.int(pvalue(object2, RET[idx], ...), diff(c(0L, idx)))

            RET <- matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
                          dimnames = dimnames(z))
            class(RET) <- c("pvalue", "matrix")
            RET
        } else {
            stop("cannot compute step-down adjusted p-values for objects of class ",
                 dQuote("NullDistribution"))
        }
    }
)

setMethod("joint",
    signature = list("MaxTypeIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, stepdown, ...) {
        if (!stepdown) {
            callNextMethod(object1, object2, stepdown, ...)
        } else {
            ## free step-down based on multivariate normality
            switch(object1@alternative,
                "less" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z)                    # smallest z first
                    pq <- length(z)
                    upper <- rep.int(Inf, pq)
                    lower <- z[o]
                },
                "greater" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z, decreasing = TRUE) # largest z first
                    pq <- length(z)
                    upper <- z[o]
                    lower <- rep.int(-Inf, pq)
                },
                "two.sided" = {
                    z <- abs(statistic(object1, type = "standardized"))
                    o <- order(z, decreasing = TRUE) # abs. largest z first
                    pq <- length(z)
                    upper <- z[o]
                    lower <- -upper
                }
            )
            Rho <- cov2cor(covariance(object1))
            RET <- numeric(pq)
            oo <- o
            for (i in 1:pq) {
                RET[i] <- pmvn(lower = lower[i], upper = upper[i],
                               mean = rep.int(0, length(oo)), corr = Rho,
                               conf.int = FALSE)
                j <- rank(oo)[1] # reindexing needed in each step
                Rho <- Rho[-j, -j, drop = FALSE]
                oo <- oo[-1]
            }
            RET <- cummax(1 - RET) # enforce monotonicity

            RET <- matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
                          dimnames = dimnames(z))
            class(RET) <- c("pvalue", "matrix")
            RET
        }
    }
)

setMethod("joint",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, stepdown, ...) {
        if (!stepdown) {
            RET <- callNextMethod(object1, object2, stepdown, ...)
        } else {
            ## free step-down based on the resampling distribution
            ## (Westfall and Young, 1993, p. 66-67, Algorithm 2.8)
            ## using standardized statistics instead of p-values
            switch(object1@alternative,
                "less" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z, decreasing = TRUE) # largest z first
                    RET <- support(object2, raw = TRUE)
                    RET <- rowMeans(colCummins(RET[o, ]) %LE% z[o])
                },
                "greater" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z)                    # smallest z first
                    RET <- support(object2, raw = TRUE)
                    RET <- rowMeans(colCummaxs(RET[o, ]) %GE% z[o])
                },
                "two.sided" = {
                    z <- abs(statistic(object1, type = "standardized"))
                    o <- order(z)                    # abs. smallest z first
                    RET <- abs(support(object2, raw = TRUE))
                    RET <- rowMeans(colCummaxs(RET[o, ]) %GE% z[o])
                }
            )
            RET <- rev(cummax(rev(RET))) # enforce monotonicity

            RET <- matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
                          dimnames = dimnames(z))
            class(RET) <- c("pvalue", "matrix")
        }

        attr(RET, "nresample") <- object2@nresample
        RET
    }
)

setMethod("joint",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, stepdown, ...) {
        callGeneric(object1@statistic, object1@distribution, stepdown, ...)
    }
)


### marginal distribution-based max-T multiple testing procedures
### (Westfall and Wolfinger, 1997; Westfall and Troendle, 2008)
setGeneric("marginal",
    function(object1, object2, ...) {
        standardGeneric("marginal")
    }
)

setMethod("marginal",
    signature = list("MaxTypeIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, stepdown, bonferroni, ...) {
        ## unadjusted p-values
        z <- statistic(object1, type = "standardized")
        RET <- switch(object1@alternative,
                   "less"      = pnorm(z),
                   "greater"   = 1 - pnorm(z),
                   "two.sided" = 2 * pmin.int(pnorm(z), 1 - pnorm(z))
               )

        ## adjustment
        RET <- if (!stepdown) {
                   if (bonferroni) # Bonferroni
                       pmin.int(1, length(RET) * RET)
                   else            # Sidak
                       1 - (1 - RET)^length(RET)
               } else {
                   n <- length(RET)
                   o <- order(RET)
                   if (bonferroni) # Bonferroni-Holm
                       pmin.int(1, cummax((n - seq_len(n) + 1L) * RET[o])[order(o)])
                   else            # Sidak-Holm
                       cummax(1 - (1 - RET[o])^(n - seq_len(n) + 1L))[order(o)]
               }

        RET <- matrix(RET, nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- c("pvalue", "matrix")
        RET
    }
)

setMethod("marginal",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, stepdown, bonferroni, ...) {
        switch(object1@alternative,
            "less" = {
                z <- -statistic(object1, type = "standardized")
                RET <- -support(object2, raw = TRUE)
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                RET <- support(object2, raw = TRUE)
            },
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                RET <- abs(support(object2, raw = TRUE))
            }
        )
        a <- attributes(z) # get dim and dimnames
        pq <- length(z)
        if (stepdown) {
            o <- order(z, decreasing = TRUE) # largest z first
            z <- z[o]
            RET <- RET[o, ]
        }

        ## marginal resampling distributions (assuming large z rejects H_0)
        RET <- lapply(1:pq, function(i) {
            Z_i <- RET[i, ]
            z_i <- unique(Z_i)
            z_i <- z_i[z_i %GE% min(z)] # optimization? '%GE%' is expensive...
            p_i <- vapply(z_i, function(z_i) mean(Z_i %GE% z_i), NA_real_)
            list(z = z_i, p = p_i)
        })
        ## single-step and free step-down based on the marginal resampling
        ## distributions (Westfall and Wolfinger, 1997) using standardized
        ## statistics instead of p-values
        RET <- vapply(1:pq, function(i) {
            p <- vapply(if (stepdown) i:pq else 1:pq, function(j) {
                RET_j <- RET[[j]]
                max(RET_j$p[RET_j$z %GE% z[i]], 0) # below eq. 2
            }, NA_real_)
            if (bonferroni) # Bonferroni(-Holm)
                pmin.int(1, sum(p)) # eq. 4
            else            # Sidak(-Holm)
                1 - prod(1 - p)     # eq. 2
        }, NA_real_)
        if (stepdown)
            RET <- cummax(RET)[order(o)] # enforce monotonicity

        attributes(RET) <- a # set dim and dimnames
        class(RET) <- c("pvalue", "matrix")
        attr(RET, "nresample") <- object2@nresample
        RET
    }
)

setMethod("marginal",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, stepdown, bonferroni, ...) {
        callGeneric(object1@statistic, object1@distribution, stepdown, bonferroni, ...)
    }
)


### unadjusted p-values
setGeneric("unadjusted",
    function(object1, object2, ...) {
        standardGeneric("unadjusted")
    }
)

setMethod("unadjusted",
    signature = list("MaxTypeIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, ...) {
        z <- statistic(object1, type = "standardized")
        RET <- switch(object1@alternative,
                   "less"      = pnorm(z),
                   "greater"   = 1 - pnorm(z),
                   "two.sided" = 2 * pmin.int(pnorm(z), 1 - pnorm(z))
               )

        RET <- matrix(RET, nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- c("pvalue", "matrix")
        RET
    }
)

setMethod("unadjusted",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, ...) {
        ## standardized observed and permuted test statistics
        switch(object1@alternative,
            "less" = {
                z <- statistic(object1, type = "standardized")
                RET <- support(object2, raw = TRUE)
                RET <- rowMeans(RET %LE% as.vector(z))
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                RET <- support(object2, raw = TRUE)
                RET <- rowMeans(RET %GE% as.vector(z))
            },
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                RET <- abs(support(object2, raw = TRUE))
                RET <- rowMeans(RET %GE% as.vector(z))
            }
        )

        RET <- matrix(RET, nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- c("pvalue", "matrix")
        attr(RET, "nresample") <- object2@nresample
        RET
    }
)

setMethod("unadjusted",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, ...) {
        callGeneric(object1@statistic, object1@distribution, ...)
    }
)


### compute p-values under subset pivotality (Westfall, 1997)
npmcp <- function(object) {

    ## extract from object
    y <- object@statistic@y[[1]]
    x <- object@statistic@x[[1]]
    ytrafo <- object@statistic@ytrafo
    alternative <- object@statistic@alternative

    ## <FIXME> it is currently hard to ask a distribution object
    ## for its type (and arguments). It's a design bug.
    distribution <- object@call$distribution
    ## </FIXME>
    z <- switch(alternative,
             "less"      = statistic(object, type = "standardized"),
             "greater"   = -statistic(object, type = "standardized"),
             "two.sided" = -abs(statistic(object, type = "standardized"))
         )

    ## get contrast matrix from xtrans
    C <- attr(object@statistic@xtrans, "contrast")
    stopifnot(inherits(C, "matrix"))

    ## order test statistics, most "extreme" one comes first
    Corder <- C[order(z), , drop = FALSE]

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

    RET <- vapply(ms, function(sub) # for every list of allowed subsets
        max(vapply(sub, foo, NA_real_)), NA_real_) # for every subset

    for (i in 2:length(RET))
        RET[i] <- max(RET[i - 1], RET[i]) # enforce monotonicity

    matrix(RET[rank(z)], dimnames = dimnames(z))
}
