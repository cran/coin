
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
    ### temporarily define `na.pass' as na.action
    opt <- options()
    on.exit(options(opt))
    options(na.action = na.pass)
    ### construct design matrix _without_ intercept
    mm <- model.matrix(~ x - 1)
    colnames(mm) <- levels(x)
    ### remove unused levels
    mm <- mm[,colSums(mm, na.rm = TRUE) > 0,drop = FALSE]
    ### the two-sample situations
    if (ncol(mm) == 2) mm <- mm[,-2,drop = FALSE]
    return(mm)
}

### transformation function
trafo <- function(data, numeric_trafo = id_trafo, factor_trafo = f_trafo, 
                 surv_trafo = logrank_trafo, var_trafo = NULL, block = NULL) {

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

    if (!is.null(var_trafo)) {
        if (!is.list(var_trafo)) stop(sQuote("var_trafo"), " is not a list")
        if (!all(names(var_trafo) %in% names(data)))
            stop("variable(s) ", 
                 names(var_trafo)[!(names(var_trafo) %in% names(data))], 
                 " not found in ", sQuote("var_trafo"))
    }

    ### compute transformations for each variable
    tr <- vector(mode = "list", length = length(data))
    names(tr) <- names(data)
    for (nm in names(data)) {
        x <- data[[nm]]
        if (nm %in% names(var_trafo)) {
            tr[[nm]] <- as.matrix(var_trafo[[nm]](x))
            next()
        }
        if (class(x)[1] == "AsIs") {
            if (length(class(x)) == 1) {
                x <- as.numeric(x)
            } else {
                class(x) <- class(x)[-1]
            }
        }
        if (is.factor(x) || is.logical(x)) {
            tr[[nm]] <- as.matrix(factor_trafo(x))
            next()
        }
        if (inherits(x, "Surv")) {
            tr[[nm]] <- as.matrix(surv_trafo(x))
            next()
        }
        if (is.numeric(x)) {
            tr[[nm]] <- as.matrix(numeric_trafo(x))
            next()
        }
        if (is.null(tr[[nm]]))
            stop("data class ", class(x), " is not supported")
    }

    ### set up a matrix of transformations
    ### when more than one factor is in play, factor names 
    ### _and_ colnames of the corresponding rows are combined by `.'
    RET <- c()
    assignvar <- c()
    cn <- c()
    for (i in 1:length(tr)) {
        if (nrow(tr[[i]]) != nrow(data)) 
            stop("Transformation of variable ", names(tr)[i], 
                 " are not of length / nrow", nrow(data))
        RET <- cbind(RET, tr[[i]])
        if (is.null(colnames(tr[[i]]))) {
            cn <- c(cn, rep("", ncol(tr[[i]])))	
        } else {
            cn <- c(cn, paste(ifelse(length(tr) > 1, ".", ""), 
                              colnames(tr[[i]]), sep = ""))
        }
        assignvar <- c(assignvar, rep(i, ncol(tr[[i]])))
    }
    attr(RET, "assign") <- assignvar
    if (length(tr) > 1) {
        colnames(RET) <- paste(rep(names(tr), tabulate(assignvar)), 
                               cn, sep = "")
    } else {
        colnames(RET) <- cn
    }
    return(RET)
}
