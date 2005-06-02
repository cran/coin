
### Streitberg-Roehmel algorithm for independent two samples
SR_shift_2sample <- function(object, fact = NULL) {

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
        stop("Argument ", sQuote("object"), " is not of class ",
             sQuote("ScalarIndependenceTestStatistic"))

    if (!is_2sample(object)) 
        stop(sQuote("object"), 
             " does not represent an independent two-sample problem")

    if (!(max(abs(object@weights - 1.0)) < sqrt(.Machine$double.eps)))
        stop("cannot compute exact distribution with non-unity weights")

    ### in case we can't map the scores into integers, use another algorithm
    if (!any(is_integer(object@ytrans[,1])))
        return(vdW_split_up_2sample(object))

    RET <- new("ExactNullDistribution")

    T <- 0
    Prob <- 1
    for (lev in levels(object@block)) {

        thisblock <- (object@block == lev)

        ### compute distribution of scores in this block
        scores <- object@ytrans[thisblock, 1]
        m <- sum(object@xtrans[thisblock, 1] == 1)

        if (m == 0) next;
        if (m == length(scores))
            dens <- list(T = sum(scores), Prob = 1)
        if (m < length(scores))
            dens <- cSR_shift_2sample(scores, m, fact = fact)

        ### update distribution of statistic over all blocks
        T <- as.vector(outer(dens$T, T, "+"))
        Prob <- drop(kronecker(Prob, dens$Prob))

    }

    T <- (T - expectation(object)) / sqrt(variance(object))

    RET@p <- function(q) sum(Prob[T <= q])
    RET@q <- function(p) {
        indx <- which(cumsum(Prob) < p)
        if (length(indx) == 0) indx <- 0
        T[max(indx) + 1]
    }
    RET@d <- function(x) Prob[T == x]
    RET@pvalue <- function(q) {
        switch(object@alternative, 
            "less"      = sum(Prob[T <= q]),
            "greater"   = sum(Prob[T >= q]),
            "two.sided" = {
                if (q == 0) return(1)
                return(sum(Prob[T <= ifelse(q >  0, -q,  q)]) + 
                       sum(Prob[T >= ifelse(q >= 0,  q, -q)]))
            }
        )
    }
    RET@support <- function(p = 1e-5) T
    RET@name <- "exact distribution (via Streitberg-Roehmel algorithm)"
    return(RET)
}

cSR_shift_2sample <- function(scores, m, fact = NULL) {

    if (m < 1 || m == length(scores))
        stop("not a two sample problem")
    n <- length(scores)
    ones <- rep(1, n)

    ### search for equivalent integer scores with sum(scores) minimal
    if (is.null(fact)) {
        fact <- c(1, 2, 10, 100, 1000)
        f <- is_integer(scores, fact = fact)
        if (!any(f))
            stop("cannot compute exact distribution with real valued scores")
        fact <- min(fact[f])
    }

    scores <- scores * fact
    add <- min(scores - 1)
    scores <- scores - add
    m_b <- sum(sort(scores)[(n + 1 - m):n])

    Prob <- .Call("R_cpermdist2", 
                  score_a = as.integer(ones),  
                  score_b = as.integer(scores),
                  m_a = as.integer(m),  
                  m_b = as.integer(m_b),
                  retProb = as.logical(TRUE), 
                  PACKAGE = "coin")

    T <- which(Prob != 0)
    Prob <- Prob[T]
    T <- (T + add*m)/fact
    return(list(T = T, Prob = Prob))
}

### van de Wiel split-up algorithm for independent two samples
vdW_split_up_2sample <- function(object) {

    ### <FIXME> on.exit(ex <- .C("FreeW", PACKAGE = "coin")) </FIXME>

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
        stop("Argument ", sQuote("object"), " is not of class ",
             sQuote("ScalarIndependenceTestStatistic"))

    if (nlevels(object@block) != 1)
        stop("cannot compute exact p-values with blocks")

    ### 2 groups as  variable
    groups <- ncol(object@xtrans) == 1 && all(object@xtrans[,1] %in% c(0, 1))
 
    RET <- new("ExactNullDistribution")

    scores <- object@ytrans[,1]
    n <- length(scores)
    storage.mode(scores) <- "double"
    m <- sum(object@xtrans)
    storage.mode(m) <- "integer"

    RET@p <- function(q) {
        obs <- q * sqrt(variance(object)) + expectation(object)
        .Call("R_split_up_2sample", scores, m, obs, PACKAGE = "coin")
    }

    RET@q <- function(p) {
        f <- function(x) RET@p(x) - p
        if (p <= 0.5)
            uniroot(f, interval = c(-10, 0))$root
        else
            uniroot(f, interval = c(0, 10))$root
    }
    RET@d <- function(x) NA

    RET@pvalue <- function(q) {
        switch(object@alternative, 
            "less"      = RET@p(q),
            "greater"   = 1 - RET@p(q - 1e-10),
            "two.sided" = {
                if (q == 0) return(1)
                if (q > 0)
                    return(RET@p(-q) + (1 - RET@p(q - 1e-10)))
                else
                    return(RET@p(q) + (1 - RET@p(- q - 1e-10)))
            }
        )
    }
    RET@support <- function(p = 1e-5) NA
    RET@name <- "exact distribution (via van de Wiel split-up algorithm)"
    return(RET)
}
