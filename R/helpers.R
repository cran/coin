
LinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_LinearStatistic", x, y, weights, 
          PACKAGE = "coin")
}

ExpectCovarInfluence <- function(y, weights) {
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_ExpectCovarInfluence", y, weights, 
          PACKAGE = "coin")
}

ExpectCovarLinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    expcovinf <- ExpectCovarInfluence(y, weights)
    .Call("R_ExpectCovarLinearStatistic", x, y, weights, expcovinf,
          PACKAGE = "coin")
}

### copied from package MASS
MPinv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X))) 
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) 
        RET <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        RET <- array(0, dim(X)[2:1])
    else RET <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
    return(list(MPinv = RET, rank = sum(Positive)))
}

copyslots <- function(source, target) {
    slots <- names(getSlots(class(source)))
    slots <- slots[(slots %in% names(getSlots(class(target))))]
    if (length(slots) == 0) 
        stop("no common slots to copy to")
    for (s in slots)
        eval(parse(text = paste("target@", s, " <- source@", s)))
    return(target)
}

formula2data <- function(formula, data, subset, ...) {

    dat <- ModelEnvFormula(formula = formula, data = data,
                           subset = subset, ...)

    ### rhs of formula
    if (has(dat, "input"))
        x <- dat@get("input")
    else 
        stop("missing right hand side of formula")

    ### ~ x + y is allowed
    if (has(dat, "response"))
        y <- dat@get("response")
    else {
        if (ncol(x) == 2) {
            y <- x[2]
            x <- x[1]
        } else 
        stop("missing left hand side of formula")
    }

    ### y ~ x | block or ~ y + x | block
    if (has(dat, "blocks")) {
        block <- dat@get("blocks")
        attr(block[[1]], "blockname") <- colnames(block)
    } else 
        block <- NULL

    ### Surv(y, event) ~ x or Surv(y, event) ~ x | block
    if (has(dat, "censored")) {
        event <- dat@get("censored")
        y <- data.frame(y = Surv(y[[1]], event[[1]]))    
    } 

    RET <- list(x = x, y = y, block = block, bl = block[[1]])

    ### <FIXME>
    if (any(sapply(RET, function(x) {
        if (is.null(x)) return(FALSE)
        if (is.list(x)) 
            return(any(sapply(x, is.na)))
        else
            return(any(is.na(x)))
    })))
        stop("NA handling currently not implemented")
    ### </FIXME>

    return(RET)
}

formula2weights <- function(weights, data, subset, ...) {

    if (is.null(weights)) return(NULL)

    if (class(weights) != "formula") 
        stop(sQuote("weights"), " is not a formula object")

    dat <- ModelEnvFormula(formula = weights, data = data,
                           subset = subset, ...)

    if (has(dat, "input"))
        weights <- dat@get("input")
    else 
        stop("missing right hand side of formula ", sQuote("weights"))

    if (length(weights) > 1 || has(dat, "response"))
        stop("only one right hand side variable allowed in ", sQuote("weights"))

    return(weights[[1]])
}

setscores <- function(x, scores) {

    if (is.null(scores)) return(x)

    if (!is.list(scores) || is.null(names(scores)))
       stop(sQuote("scores"), " is not a named list")

    varnames <- names(scores)

    missing <- varnames[!varnames %in% c(colnames(x@x), colnames(x@y))]
    if (length(missing) > 0)
        stop("Variable(s)", paste(missing, sep = ", "), 
             " not found in ", sQuote("x"))


    for (var in varnames) {
        if (!is.null(x@x[[var]])) {

            if (!is.factor(x@x[[var]]))
                stop(var, " is not a factor")

            if (nlevels(x@x[[var]]) != length(scores[[var]]))
                stop("scores for variable ", var, " don't match")
            
            x@x[[var]] <- ordered(x@x[[var]], levels = levels(x@x[[var]]))
            attr(x@x[[var]], "scores") <- scores[[var]]
        }
        if (!is.null(x@y[[var]])) {

            if (!is.factor(x@y[[var]]))
                stop(var, " is not a factor")

            if (nlevels(x@y[[var]]) != length(scores[[var]]))
                stop("scores for variable ", var, " don't match")
            
            x@y[[var]] <- ordered(x@y[[var]], levels = levels(x@y[[var]]))
            attr(x@y[[var]], "scores") <- scores[[var]]
        }
    }
    return(x)
}

check_trafo <- function(tx, ty) {

    if (!is.matrix(tx))
        stop(sQuote("xtrafo"), " does not return a matrix")
    if (!is.matrix(ty))
        stop(sQuote("ytrafo"), " does not return a matrix")
    if (nrow(tx) != nrow(ty))
        stop("Dimensions of returns of ", sQuote("xtrafo"), " and ",
             sQuote("ytrafo"), " don't match")
    storage.mode(tx) <- "double"
    storage.mode(ty) <- "double"
    list(xtrafo = tx, ytrafo = ty)
}

table2df <- function(x) {
    if (!is.table(x))
        stop(sQuote("x"), " is not of class ", sQuote("table"))
    x <- as.data.frame(x)
    freq <- x[["Freq"]]
    x <- x[rep(1:nrow(x), freq), ,drop = FALSE]
    rownames(x) <- 1:nrow(x)
    return(x[,colnames(x) != "Freq"])
}

table2df_sym <- function(x) {
    x <- table2df(x)
    lx <- levels(x[[1]])
    if (!all(sapply(x, function(x) all(levels(x) == lx))))
        stop("table ", sQuote("x"), " does not represent a symmetry problem")
    n <- nrow(x)
    p <- ncol(x)
    y <- data.frame(groups = factor(rep(colnames(x), rep(n, p))),
                    response = factor(unlist(x), labels = lx))
    rownames(y) <- 1:(n*p)
    y
}

is_2sample <- function(object) {
    groups <- ((ncol(object@x) == 1 && is.factor(object@x[[1]])) && 
                nlevels(object@x[[1]]) == 2)
    groups <- groups && ncol(object@xtrans) == 1
    values <- (ncol(object@y) == 1 && ncol(object@ytrans) == 1)
    return(groups && values)
}

is_Ksample <- function(object) {
    groups <- (ncol(object@x) == 1 && is.factor(object@x[[1]]))
    values <- (ncol(object@y) == 1 && ncol(object@ytrans) == 1)
    return(groups && values)
}

is_numeric_y <- function(object) {
    is.numeric(object@y[[1]])
}

is_censored_y <- function(object) {
    ncol(object@y) == 1 && class(object@y[[1]]) == "Surv"
}

is_corr <- function(object) {
    (is.numeric(object@x[[1]]) && is.numeric(object@y[[1]])) &&
     (ncol(object@xtrans) == 1 && ncol(object@ytrans) == 1)
}

is_contingency <- function(object) {
    x <- object@x
    groups <- (ncol(x) == 1 && is.factor(x[[1]]))
    y <- object@y
    values <- (ncol(y) == 1 && is.factor(y[[1]]))
    trans <- all(rowSums(object@xtrans) %in% c(0,1)) && 
             all(rowSums(object@ytrans) %in% c(0, 1))
    return((groups && values) && trans)       
}

is_ordered <- function(object) {
    x <- object@x
    y <- object@y
    (is_Ksample(object) || is_contingency(object)) && 
    (is.ordered(x[[1]]) || is.ordered(y[[1]]))
}

is_completeblock <- function(object) {
    x <- object@x
    all(table(object@x[[1]], object@block) == 1)
}

is_scalar <- function(object) {
    ncol(object@xtrans) == 1 && ncol(object@ytrans) == 1
}

is_ordered_x <- function(object) {
    all(sapply(object@x, function(x) is.numeric(x) || is.ordered(x)))
}

is_integer <- function(x, fact = c(1, 2, 10, 100, 1000))
    sapply(fact, function(f) max(abs(round(x * f) - (x * f))) < 
           sqrt(.Machine$double.eps))

is_censored <- function(object) {
    ncol(object@y) == 1 && class(object@y[[1]]) == "Surv"
}

isequal <- function(a, b) {
    attributes(a) <- NULL
    if (!identical(all.equal(a, b), TRUE)) {
        print(a, digits = 10)
        print(b, digits = 10)
        return(FALSE)
    } else {
        return(TRUE)
    }
}
