### Streitberg-Roehmel algorithm for two independent samples
SR_shift_2sample <- function(object, fact) {
    teststat <-
        if (inherits(object, "ScalarIndependenceTestStatistic"))
            "scalar"
        else if (inherits(object, "QuadTypeIndependenceTestStatistic"))
            "quadratic"
        else
            stop(sQuote("object"), " is not of class ",
                 dQuote("ScalarIndependenceTestStatistic"), " or ",
                 dQuote("QuadTypeIndependenceTestStatistic"))

    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent an independent two-sample problem")

    ytrans <- object@ytrans[, 1L]
    xtrans <- object@xtrans[, 1L]
    block <- object@block

    ## expand observations if weights are non-unity
    if (!is_unity(object@weights)) {
        idx <- rep.int(seq_along(object@weights), object@weights)
        ytrans <- ytrans[idx]
        xtrans <- xtrans[idx]
        block <- block[idx]
    }

    lev <- levels(block)
    nb <- nlevels(block)

    ## first block
    block_1 <- (block == lev[1L])
    scores <- ytrans[block_1]
    m <- sum(xtrans[block_1] == 1L)

    ## compute T and density
    if (m == 0L)
        dens <- list(T = 0, Prob = 1)
    else if (m == length(scores))
        dens <- list(T = sum(scores), Prob = 1)
    else if (m < length(scores))
        dens <- cSR_shift_2sample(scores, m, fact = fact)
    else
        stop("cannot compute exact distribution")

    T <- dens$T
    Prob <- dens$Prob

    ## remaining blocks
    if (nb > 1) {
        for (i in seq_len(nb)[-1L]) {
            block_i <- (block == lev[i])
            scores <- ytrans[block_i]
            m <- sum(xtrans[block_i] == 1L)

            ## compute T and density
            if (m == 0L)
                next
            else if (m == length(scores))
                dens <- list(T = sum(scores), Prob = 1)
            else if (m < length(scores))
                dens <- cSR_shift_2sample(scores, m, fact = fact)
            else
                stop("cannot compute exact distribution")

            ## update T
            T <- .Call(R_outersum, dens$T, T)

            ## make sure T is ordered and distinct
            n <- length(T)
            o <- order(T)
            T <- T[o]
            idx <- c(which(T[-1L] %NE% T[-n]), n)
            T <- T[idx]

            ### update density (use C_kronecker from libcoin)
            Prob <- .Call(R_kronecker, dens$Prob, Prob)
            Prob <- vapply(split(Prob[o],
                                 rep.int(seq_along(idx), diff(c(0L, idx)))),
                           sum, NA_real_, USE.NAMES = FALSE)
        }
    }

    if (teststat == "scalar")
        T <- (T - expectation(object)) / sqrt(variance(object))
    else {
        T <- (T - expectation(object))^2 / variance(object)
        ## make sure T is ordered and distinct
        n <- length(T)
        o <- order(T)
        T <- T[o]
        idx <- c(which(T[-1L] %NE% T[-n]), n)
        T <- T[idx]
        ## compute density
        Prob <- vapply(split(Prob[o],
                             rep.int(seq_along(idx), diff(c(0L, idx)))),
                       sum, NA_real_, USE.NAMES = FALSE)
    }

    p_fun <- function(q) {
        sum(Prob[T %LE% q])
    }
    q_fun <- function(p) {
        idx <- which(cumsum(Prob) %LT% p)
        if (length(idx) == 0L)
            T[1L]
        else if (length(idx) == length(Prob))
            T[max(idx)]
        else
            T[max(idx) + 1L]
    }
    d_fun <- function(x) {
        eq <- T %EQ% x
        if (any(eq))
            Prob[eq]
        else
            0
    }
    pvalue_fun <- function(q) {
        if (teststat == "scalar")
            switch(object@alternative,
                "less"      = sum(Prob[T %LE% q]),
                "greater"   = sum(Prob[T %GE% q]),
                "two.sided" = {
                    if (q == 0)
                        1L
                    else
                        sum(Prob[T %LE% if (q < 0) q else -q]) +
                          sum(Prob[T %GE% if (q > 0) q else -q])
                }
            )
        else {
            if (q == 0)
                1L
            else
                sum(Prob[T %GE% q])
        }
    }
    midpvalue_fun <- function(q, z) {
        pvalue_fun(q) - z *
          if (teststat == "scalar" && object@alternative == "two.sided")
              d_fun(-q) + d_fun(q) # both tails
          else
              d_fun(q)
    }

    p <- function(q) {
        ## NOTE: 'vapply' provide names
        vapply(q, p_fun, NA_real_)
    }
    q <- function(p) {
        ## NOTE: 'vapply' provide names
        vapply(p, q_fun, NA_real_)
    }
    d <- function(x) {
        ## NOTE: 'vapply' provide names
        vapply(x, d_fun, NA_real_)
    }
    pvalue <- function(q) {
        if (length(q) < 2L)
            pvalue_fun(q)
        else
            vapply(q, pvalue_fun, NA_real_)
    }
    midpvalue <- function(q) {
        if (length(q) < 2L)
            midpvalue_fun(q, z = 0.5)
        else
            vapply(q, midpvalue_fun, NA_real_, z = 0.5)
    }
    pvalueinterval <- function(q) {
        if (length(q) < 2L)
            midpvalue_fun(q, z = c("p_0" = 1, "p_1" = 0))
        else
            vapply(q, midpvalue_fun, c(NA_real_, NA_real_),
                   z = c("p_0" = 1, "p_1" = 0))
    }
    support <- function() T
    size <- function(alpha, type) {
        pv_fun <- if (type == "mid-p-value") midpvalue else pvalue
        spt <- support()
        vapply(alpha, function(a) sum(d(spt[pv_fun(spt) %LE% a])), NA_real_)
    }

    new("ExactNullDistribution",
        p = p,
        q = q,
        d = d,
        pvalue = pvalue,
        midpvalue = midpvalue,
        pvalueinterval = pvalueinterval,
        size = size,
        support = support,
        name = paste0("Exact Distribution for Independent Two-Sample Tests",
                      " (Streitberg-Roehmel Shift Algorithm)"))
}

cSR_shift_2sample <- function(scores, m, fact) {
    n <- length(scores)
    if (m < 1L || m == n)
        stop("not a two sample problem")

    ## integer scores with sum(scores) minimal
    scores <- round(scores * fact)
    storage.mode(scores) <- "integer"
    add <- min(scores - 1L)
    scores <- scores - add
    storage.mode(m) <- "integer"
    m_b <- sum(sort(scores)[(n + 1L - m):n])

    Prob <- .Call(R_cpermdist2,
                  score_a = rep.int(1L, n), score_b = scores,
                  m_a = m, m_b = m_b, retProb = TRUE)
    T <- which(Prob != 0)

    list(T = (T + add * m) / fact, Prob = Prob[T])
}


### Streitberg-Roehmel algorithm for two paired samples
SR_shift_1sample <- function(object, fact) {
    teststat <-
        if (inherits(object, "ScalarIndependenceTestStatistic"))
            "scalar"
        else if (inherits(object, "QuadTypeIndependenceTestStatistic"))
            "quadratic"
        else
            stop(sQuote("object"), " is not of class ",
                 dQuote("ScalarIndependenceTestStatistic"), " or ",
                 dQuote("QuadTypeIndependenceTestStatistic"))

    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent an independent two-sample problem")

    scores <- object@ytrans[, 1L]
    if (any(scores < 0))
        stop("cannot compute exact distribution with negative scores")
    block <- object@block

    ## expand observations if weights are non-unity
    if (!is_unity(object@weights)) {
        idx <- rep.int(seq_along(object@weights), object@weights)
        scores <- scores[idx]
        block <- block[idx]
    }

    ##  table(object@block, scores == 0) checken
    scores <- vapply(unique(object@block), function(i) {
        s <- round(scores * fact)[object@block == i]
        s[s != 0] # remove zeros
    }, NA_real_)
    storage.mode(scores) <- "integer"
    Prob <- .Call(R_cpermdist1, scores)
    T <- which(Prob != 0)
    Prob <- Prob[T]
    ## 0 is possible
    T <- (T - 1) / fact

    if (teststat == "scalar")
        T <- (T - expectation(object)) / sqrt(variance(object))
    else {
        T <- (T - expectation(object))^2 / variance(object)
        ## make sure T is ordered and distinct
        n <- length(T)
        o <- order(T)
        T <- T[o]
        idx <- c(which(T[-1L] %NE% T[-n]), n)
        T <- T[idx]
        ## compute density
        Prob <- vapply(split(Prob[o],
                             rep.int(seq_along(idx), diff(c(0L, idx)))),
                       sum, NA_real_, USE.NAMES = FALSE)
    }

    p_fun <- function(q) {
        sum(Prob[T %LE% q])
    }
    q_fun <- function(p) {
        idx <- which(cumsum(Prob) %LT% p)
        if (length(idx) == 0L)
            T[1L]
        else if (length(idx) == length(Prob))
            T[max(idx)]
        else
            T[max(idx) + 1L]
    }
    d_fun <- function(x) {
        eq <- T %EQ% x
        if (any(eq))
            Prob[eq]
        else
            0
    }
    pvalue_fun <- function(q) {
        if (teststat == "scalar")
            switch(object@alternative,
                "less"      = sum(Prob[T %LE% q]),
                "greater"   = sum(Prob[T %GE% q]),
                "two.sided" = {
                    if (q == 0)
                        1L
                    else
                        sum(Prob[T %LE% if (q < 0) q else -q]) +
                          sum(Prob[T %GE% if (q > 0) q else -q])
                }
            )
        else {
            if (q == 0)
                1L
            else
                sum(Prob[T %GE% q])
        }
    }
    midpvalue_fun <- function(q, z) {
        pvalue_fun(q) - z *
          if (teststat == "scalar" && object@alternative == "two.sided")
              d_fun(-q) + d_fun(q) # both tails
          else
              d_fun(q)
    }

    p <- function(q) {
        ## NOTE: 'vapply' provide names
        vapply(q, p_fun, NA_real_)
    }
    q <- function(p) {
        ## NOTE: 'vapply' provide names
        vapply(p, q_fun, NA_real_)
    }
    d <- function(x) {
        ## NOTE: 'vapply' provide names
        vapply(x, d_fun, NA_real_)
    }
    pvalue <- function(q) {
        if (length(q) < 2L)
            pvalue_fun(q)
        else
            vapply(q, pvalue_fun, NA_real_)
    }
    midpvalue <- function(q) {
        if (length(q) < 2L)
            midpvalue_fun(q, z = 0.5)
        else
            vapply(q, midpvalue_fun, NA_real_, z = 0.5)
    }
    pvalueinterval <- function(q) {
        if (length(q) < 2L)
            midpvalue_fun(q, z = c("p_0" = 1, "p_1" = 0))
        else
            vapply(q, midpvalue_fun, c(NA_real_, NA_real_),
                   z = c("p_0" = 1, "p_1" = 0))
    }
    support <- function() T
    size <- function(alpha, type) {
        pv_fun <- if (type == "mid-p-value") midpvalue else pvalue
        spt <- support()
        vapply(alpha, function(a) sum(d(spt[pv_fun(spt) %LE% a])), NA_real_)
    }

    new("ExactNullDistribution",
        p = p,
        q = q,
        d = d,
        pvalue = pvalue,
        midpvalue = midpvalue,
        pvalueinterval = pvalueinterval,
        size = size,
        support = support,
        name = paste0("Exact Distribution for Dependent Two-Sample Tests",
                      " (Streitberg-Roehmel Shift Algorithm)"))
}


### van de Wiel split-up algorithm for two independent samples
vdW_split_up_2sample <- function(object) {
    ## <FIXME> on.exit(ex <- .C("FreeW", PACKAGE = "coin")) </FIXME>

    if (!inherits(object, "ScalarIndependenceTestStatistic"))
        stop(sQuote("object"), " is not of class ",
             dQuote("ScalarIndependenceTestStatistic"))

    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent an independent two-sample problem")

    if (nlevels(object@block) != 1L)
        stop("cannot compute exact p-values with blocks")

    scores <- object@ytrans[, 1L]
    xtrans <- object@xtrans[, 1L]

    ## expand observations if weights are non-unity
    if (!is_unity(object@weights)) {
        idx <- rep.int(seq_along(object@weights), object@weights)
        scores <- scores[idx]
        xtrans <- xtrans[idx]
    }

    storage.mode(scores) <- "double"
    m <- sum(xtrans)
    storage.mode(m) <- "integer"

    p_fun <- function(q) {
        T <- q * sqrt(variance(object)) + expectation(object)
        .Call(R_split_up_2sample, scores, m, T, sqrt_eps)
    }
    q_fun <- function(p) {
        f <- function(x) p_fun(x) - p
        rr <- if (p <= 0.5)
                  uniroot(f, interval = c(-10, 1), tol = sqrt_eps)
              else
                  uniroot(f, interval = c(-1, 10), tol = sqrt_eps)
        ## make sure quantile leads to pdf >= p
        if (rr$f.root < 0)
            rr$root <- rr$root + sqrt_eps
        ## pdf is constant here
        if (rr$estim.prec > sqrt_eps) {
            r1 <- rr$root
            d <- min(diff(sort(scores[!duplicated(scores)]))) /
                   sqrt(variance(object))
            while (d > sqrt_eps) {
                if (f(r1 - d) >= 0)
                    r1 <- r1 - d
                else
                    d <- d / 2
            }
            rr$root <- r1
        }
        rr$root
    }
    pvalue_fun <- function(q) {
        switch(object@alternative,
            "less"      = p_fun(q),
            "greater"   = 1 - p_fun(q - 10 * sqrt_eps),
            "two.sided" = {
                if (q == 0)
                    1L
                else if (q > 0)
                    p_fun(-q) + (1 - p_fun(q - 10 * sqrt_eps))
                else
                    p_fun(q) + (1 - p_fun(-q - 10 * sqrt_eps))
            }
        )
    }

    p <- function(q) {
        ## NOTE: 'vapply' provide names
        vapply(q, p_fun, NA_real_)
    }
    q <- function(p) {
        ## NOTE: 'vapply' provide names
        vapply(p, q_fun, NA_real_)
    }
    pvalue <- function(q) {
        if (length(q) < 2L)
            pvalue_fun(q)
        else
            vapply(q, pvalue_fun, NA_real_)
    }

    new("ExactNullDistribution",
        p = p,
        q = q,
        d = function(x) NA,
        pvalue = pvalue,
        midpvalue = function(q) NA,
        pvalueinterval = function(q) NA,
        size = function(alpha, type) NA,
        support = function() NA,
        name = paste0("Exact Distribution for Independent Two-Sample Tests",
                      " (van de Wiel Split-Up Algorithm)"))
}
