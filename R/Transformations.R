
### compute average scores, see Hajek, Sidak, Sen (page 131ff)
average_scores <- function(s, x) {
    dup <- x[duplicated(x)]
    for (d in dup)
        s[x == d] <- mean(s[x == d])
    return(s)
} 

### identity transformation
id_trafo <- function(x) x

### Ansari-Bradley
ansari_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method, 
        "mid-ranks" = {
            r <- rank(x)
            pmin(r, length(x) - r + 1)
        },
        "average-scores" = {
            r <- rank(x, ties.method = "random")
            s <- pmin(r, length(x) - r + 1)
            average_scores(s, x)
        }
    )
    return(scores)
}

### Fligner
fligner_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method, 
        "mid-ranks" = {
            qnorm((1 + rank(abs(x))/(length(x) + 1))/2)
        },
        "average-scores" = {
            r <- rank(abs(x), ties.method = "random")
            s <- qnorm((1 + r/(length(x) + 1))/2)
            average_scores(s, x)
        }
    )
    return(scores)
}

### Normal Scores (van der Waerden)
normal_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = { 
            qnorm(rank(x)/(length(x) + 1))
        },
        "average-scores" = {
            r <- rank(x, ties.method = "random")
            s <- qnorm(r/(length(x) + 1))
            average_scores(s, x)
        }
    )
    return(scores)
}
 
### Median Scores
median_trafo <- function(x)
    as.numeric(x <= median(x))

### Conover & Salsburg (1988)
consal_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    scores <- switch(ties.method,
        "mid-ranks" = { 
            (rank(x)/(length(x) + 1))^4
        },
        "average-scores" = {
            r <- rank(x, ties.method = "random")
            s <- (r/(length(x) + 1))^4
            average_scores(s, x)
        }
    )
    return(scores)
}

### maximally selected (rank, chi^2, whatsoever) statistics
maxstat_trafo <- function(x, minprob = 0.1, maxprob = 0.9) {
    qx <- quantile(x, prob = c(minprob, maxprob), type = 1)
    ux <- unique(x)
    ux <- ux[ux < max(x) & ux > min(x)]
    cutpoints <- ux[ux > qx[1] & ux <= qx[2]]
    cm <- .Call("R_maxstattrafo", as.double(x), as.double(cutpoints),
                PACKAGE = "coin")
    colnames(cm) <- paste("x <= ", round(cutpoints, 3), sep = "")
    rownames(cm) <- 1:nrow(cm)
    cm
}

### logrank scores; with two different methods of handling
### ties
logrank_trafo <- function(x, ties.method = c("logrank", "HL")) {
    ties.method <- match.arg(ties.method)
    time <- x[,1]
    event <- x[,2]
    n <- length(time)
    ot <- order(time, event)
    rt <- rank(time, ties.method = "max")
    mt <- rank(time, ties.method = "min") - 1
    fact <- switch(ties.method, "logrank" = event / (n - mt),
                                "HL" = event/(n - rt + 1)
                  )
    event - cumsum(fact[ot])[rt]
}

### factor handling
f_trafo <- function(x) {
    mm <- model.matrix(~ x - 1)
    colnames(mm) <- levels(x)
    mm <- mm[,colSums(mm) > 0,drop = FALSE]
    ### the two-sample situations
    if (ncol(mm) == 2) mm <- mm[,-2,drop = FALSE]
    return(mm)
}

### transformation function
trafo <- function(data, numeric_trafo = id_trafo, factor_trafo = f_trafo, 
                 surv_trafo = logrank_trafo, block = NULL) {

    if (!(is.data.frame(data) || is.list(data)))
        stop(sQuote("data"), " is not a data.frame or list")

    if (!is.null(block)) {
        if (!is.factor(block) || length(block) != nrow(data))
            stop(sQuote("block"), " is not a factor with ", 
                 nrow(data), " elements")

        ### need to check dimension of matrix returned by 
        ### user supplied functions
        ret <- trafo(data, numeric_trafo, factor_trafo, surv_trafo)

        ### apply trafo to each block separately
        for (lev in levels(block)) {
            ret[block == lev, ] <- trafo(data[block == lev, ,drop = FALSE], 
                numeric_trafo, factor_trafo, surv_trafo)
        }
        return(ret)
    }

    tr <- lapply(data, function(x) {
        if (is.factor(x))
            return(factor_trafo(x))
        if (class(x) == "Surv")
            return(surv_trafo(x))
        if (is.numeric(x))
            return(numeric_trafo(x))
        stop("data class ", class(x), " is not supported")
    })

    chk <- sapply(tr, function(x) (is.matrix(x) && nrow(x) == nrow(data)) ||
                           (is.vector(x) && length(x) == nrow(data)))
    if (!all(chk))
        stop("transformations are not of length / nrow", nrow(data))

    RET <- c()
    assignvar <- c()
    for (i in 1:length(tr)) {
        RET <- cbind(RET, tr[[i]])
        p <- ifelse(is.matrix(tr[[i]]), ncol(tr[[i]]), 1)
        assignvar <- c(assignvar, rep(i, p))
    }
    attr(RET, "assign") <- assignvar
    return(RET)
}
