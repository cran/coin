asymptotic <- function(maxpts = 25000, abseps = 0.001, releps = 0) {
    function(object)
        AsymptNullDistribution(object, maxpts = maxpts, abseps = abseps,
                               releps = releps)
}

approximate <- function(nresample = 10000L, parallel = c("no", "multicore", "snow"),
                        ncpus = 1L, cl = NULL, B) {
    ## <DEPRECATED>
    if (!missing(B)) {
        warning(sQuote("B"), " is deprecated; use ", sQuote("nresample"),
                " instead")
        nresample <- B
    }
    ## </DEPRECATED>
    parallel <- match.arg(parallel)
    function(object)
        ApproxNullDistribution(object, nresample = nresample, parallel = parallel,
                               ncpus = ncpus, cl = cl)
}

exact <- function(algorithm = c("auto", "shift", "split-up"), fact = NULL) {
    algorithm <- match.arg(algorithm)
    function(object)
        ExactNullDistribution(object, algorithm = algorithm, fact = fact)
}

pmvn <- function(lower, upper, mean, corr, conf.int, ...) {
    p <- if (length(corr) > 1L)
             pmvnorm(lower = lower, upper = upper, mean = mean,
                     corr = corr, ...)
         else
             pmvnorm(lower = lower, upper = upper, mean = mean,
                     sigma = 1, ...)
    if (conf.int) {
        error <- attr(p, "error")
        ci <- c(max(0, p - error), min(p + error, 1))
        attributes(ci) <- list("conf.level" = 0.99)
        attributes(p) <- list("conf.int" = ci)
    } else
        attributes(p) <- NULL
    p
}

qmvn <- function(p, mean, corr, ...) {
    q <- if (length(corr) > 1L)
             qmvnorm(p = p, mean = mean, corr = corr,
                     tail = "both.tails", ...)$quantile
         else
             qmvnorm(p = p, mean = mean, sigma = 1,
                     tail = "both.tails", ...)$quantile
    attributes(q) <- NULL
    q
}

### copied from package MASS
MPinv <- function (X, tol = sqrt_eps)
{
    if (length(dim(X)) > 2L || !is.numeric(X))
        stop(sQuote("X"), " must be a numeric matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    RET <- if (all(Positive))
               Xsvd$v %*% (1 / Xsvd$d * t(Xsvd$u))
           else if (!any(Positive))
               array(0, dim(X)[2L:1L])
           else
               Xsvd$v[, Positive, drop = FALSE] %*%
                 ((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    list(MPinv = RET, rank = sum(Positive))
}

copyslots <- function(source, target) {
    slots <- slotNames(source)
    slots <- slots[slots %in% slotNames(target)]
    if (length(slots) == 0)
        stop("no common slots to copy to")
    eval(str2expression(paste0("target@", slots, " <- source@", slots)))
    target
}

ft <- function(name, class, formula, data = list(), subset = NULL,
               weights = NULL, ...) {

    object <- formula2data(formula, data, subset, weights = weights, ...)
    object <- new(class, x = object$x, y = object$y, block = object$block,
                  weights = object$weights)
    args <- list(...)
    args$frame <- NULL

    ## warn users of weighted rank tests
    if (name %in% ranktests && !is.null(object@weights) &&
        !is_unity(object@weights))
        warning("rank transformation doesn't take weights into account")

    do.call(name, c(object = object, args))
}

ranktests <-
    c("wilcox_test", "kruskal_test", "normal_test", "median_test",
      "savage_test", "taha_test", "klotz_test", "mood_test", "ansari_test",
      "fligner_test", "conover_test", "logrank_test", "quade_test",
      "friedman_test", "wilcoxsign_test", "spearman_test", "fisyat_test",
      "quadrant_test", "koziol_test")

formula2data <- function(formula, data, subset, weights = NULL, ...) {
    no_weights <- is.null(weights)

    dat <- ModelEnvFormula(
        formula = formula,
        data = data,
        subset = subset,
        other = if (no_weights) list() else list(weights = weights),
        na.action = na.omit,
        designMatrix = FALSE, responseMatrix = FALSE,
        ...
    )

    ## rhs of formula
    if (has(dat, "input"))
        x <- dat@get("input")
    else
        stop("missing right hand side of formula")

    ## ~ x + y is allowed
    if (has(dat, "response"))
        y <- dat@get("response")
    else {
        if (ncol(x) == 2L) {
            y <- x[2L]
            x <- x[1L]
        } else
            stop("missing left hand side of formula")
    }

    ## y ~ x | block or ~ y + x | block
    if (has(dat, "blocks")) {
        block <- dat@get("blocks")
        attr(block[[1L]], "blockname") <- colnames(block)
    } else
        block <- NULL

    list(x = x, y = y, block = block[[1L]],
         weights = if (no_weights) NULL else dat@get("weights")[[1L]])
}

setscores <- function(x, scores) {

    if (is.null(scores)) return(x)

    varnames <- names(scores)
    if (!is.list(scores) || is.null(varnames))
        stop(sQuote("scores"), " is not a named list")

    missing <- varnames[!varnames %in% c(colnames(x@x), colnames(x@y))]
    if (length(missing) > 0L)
        stop("variable(s) ", paste(missing, sep = ", "),
             " not found in ", sQuote("x"))
    ## <FIXME> Repeated reassignment to S4 objects may be expensive
    for (var in varnames) {
        if (!is.null(x@x[[var]])) {
            if (!is.factor(x@x[[var]]))
                stop(sQuote(var), " is not a factor")
            if (nlevels(x@x[[var]]) != length(scores[[var]]))
                stop("scores for variable ", sQuote(var), " don't match")
            x@x[[var]] <- ordered(x@x[[var]], levels = levels(x@x[[var]]))
            attr(x@x[[var]], "scores") <- scores[[var]]
        }
        if (!is.null(x@y[[var]])) {
            if (!is.factor(x@y[[var]]))
                stop(sQuote(var), " is not a factor")
            if (nlevels(x@y[[var]]) != length(scores[[var]]))
                stop("scores for variable ", sQuote(var), " don't match")
            x@y[[var]] <- ordered(x@y[[var]], levels = levels(x@y[[var]]))
            attr(x@y[[var]], "scores") <- scores[[var]]
        }
    }
    ## </FIXME>
    x
}

### user-supplied trafo functions may return a vector or matrices
### with NROW being equal for the x and y variables
check_trafo <- function(tx, ty) {

    if (!(is.numeric(tx) || is.logical(tx)))
        stop(sQuote("xtrafo"), " does not return a numeric or logical vector")
    if (!(is.numeric(ty) || is.logical(ty)))
        stop(sQuote("ytrafo"), " does not return a numeric or logical vector")
    if (NROW(tx) != NROW(ty))
        stop("dimensions of returns of ", sQuote("xtrafo"), " and ",
             sQuote("ytrafo"), " don't match")
    if (!is.matrix(tx)) tx <- matrix(tx, ncol = 1L)
    if (!is.matrix(ty)) ty <- matrix(ty, ncol = 1L)
    storage.mode(tx) <- "double"
    storage.mode(ty) <- "double"
    list(xtrafo = tx, ytrafo = ty)
}

table2df <- function(x) {
    if (!is.table(x))
        stop(sQuote("x"), " is not of class ", dQuote("table"))
    x <- as.data.frame(x)
    freq <- x[["Freq"]]
    x <- x[rep.int(seq_len(nrow(x)), freq), , drop = FALSE]
    rownames(x) <- NULL
    x[, colnames(x) != "Freq"]
}

table2df_sym <- function(x) {
    x <- table2df(x)
    lx <- levels(x[[1L]])
    if (!all(vapply(x, function(x) all(levels(x) == lx), NA)))
        stop("table ", sQuote("x"), " does not represent a symmetry problem")
    data.frame(conditions = factor(rep.int(seq_len(ncol(x)),
                                           rep.int(nrow(x), ncol(x))),
                                   labels = colnames(x)),
               response = factor(unlist(x, recursive = FALSE,
                                        use.names = FALSE),
                                 labels = lx))
}

table2IndependenceProblem <-
    function(x)
{
    x <- as.data.frame(x)
    if (ncol(x) == 3L)
        new("IndependenceProblem",
            x = x[1L], y = x[2L], block = NULL, weights = x[["Freq"]])
    else if (ncol(x) == 4L) {
        attr(x[[3L]], "blockname") <- colnames(x)[3L]
        new("IndependenceProblem",
            x = x[1L], y = x[2L], block = x[[3L]], weights = x[["Freq"]])
    } else
        stop(sQuote("object"), " is not a two- or three-way contingency table")
}

table2SymmetryProblem <-
    function(x)
{
    x <- as.data.frame(x)
    ## SymmetryProblem cannot handle weights, so expand manually
    x <- x[rep.int(seq_len(nrow(x)), x[["Freq"]]), -ncol(x)]
    lx <- levels(x[[1L]])
    if (!all(vapply(x, function(x) all(levels(x) == lx), NA)))
        stop("table ", sQuote("x"), " does not represent a symmetry problem")
    new("SymmetryProblem",
        x = data.frame(conditions = gl(ncol(x), nrow(x), labels = colnames(x))),
        y = data.frame(response = unlist(x, recursive = FALSE, use.names = FALSE)))
}


is_ytrafo <- function()
    any(vapply(sys.calls(), function(i)
            identical(as.character(i)[1], "ytrafo"), NA))

is_factor_y <- function(object)
    ncol(object@y) == 1L && is.factor(object@y[[1L]])

is_factor_x <- function(object)
    ncol(object@x) == 1L && is.factor(object@x[[1L]])

is_ordered_y <- function(object)
    ncol(object@y) == 1L && is.ordered(object@y[[1L]])

is_ordered_x <- function(object)
    ncol(object@x) == 1L && is.ordered(object@x[[1L]])

is_unordered_y <- function(object)
    ncol(object@y) == 1L && is.factor(object@y[[1L]]) && !is.ordered(object@y[[1L]])

is_unordered_x <- function(object)
    ncol(object@x) == 1L && is.factor(object@x[[1L]]) && !is.ordered(object@x[[1L]])

is_numeric_y <- function(object)
    ncol(object@y) == 1L && is.numeric(object@y[[1L]]) && !is.Surv(object@y[[1L]])

is_numeric_x <- function(object)
    ncol(object@x) == 1L && is.numeric(object@x[[1L]]) && !is.Surv(object@x[[1L]])

is_censored_y <- function(object)
    ncol(object@y) == 1L && is.Surv(object@y[[1L]])

is_censored_x <- function(object)
    ncol(object@x) == 1L && is.Surv(object@x[[1L]])

is_2sample <- function(object)
    ncol(object@x) == 1L && nlevels(object@x[[1L]]) == 2L

is_Ksample <- is_factor_x

is_corr <- function(object)
    (is_numeric_y(object) || is_censored_y(object)) &&
        (is_numeric_x(object) || is_censored_x(object))

is_contingency <- function(object)
    is_factor_y(object) && is_factor_x(object)

is_contingency_2x2 <- function(object)
    is_contingency(object) && nlevels(object@y[[1L]]) == 2L &&
        nlevels(object@x[[1L]]) == 2L

is_singly_ordered <- function(object) {
    ## NOTE: unordered factors with exactly 2 levels are regarded as ordered
    (is_ordered_y(object) &&
     (is_numeric_x(object) || is_censored_x(object) ||
      (is_unordered_x(object) && nlevels(object@x[[1L]]) > 2L))) ||
    (is_ordered_x(object) &&
     (is_numeric_y(object) || is_censored_y(object) ||
      (is_unordered_y(object) && nlevels(object@y[[1L]]) > 2L)))
}

is_doubly_ordered <- function(object) {
    ## NOTE: unordered factors with exactly 2 levels are regarded as ordered
    (is_ordered_y(object) &&
     (is_ordered_x(object) ||
      (is_unordered_x(object) && nlevels(object@x[[1L]]) == 2L))) ||
    (is_ordered_x(object) &&
     (is_ordered_y(object) ||
      (is_unordered_y(object) && nlevels(object@y[[1L]]) == 2L)))
}

is_ordered <- function(object)
    is_singly_ordered(object) || is_doubly_ordered(object)

is_completeblock <- function(object)
    all(table(object@x[[1L]], object@block) == 1L)

is_scalar <- function(object)
    ncol(object@xtrans) == 1L && ncol(object@ytrans) == 1L

is_integer <- function(x, fact = NULL) {
    if (is.null(fact))
        fact <- c(1, 2, 10, 100, 1000, 10000, 100000)
    f <- vapply(fact, function(f) max(abs(round(x * f) - (x * f))) < sqrt_eps, NA)
    if (RET <- any(f))
        attr(RET, "fact") <- min(fact[f])
    RET
}

is_monotone <- function(x)
    all(x == cummax(x)) || all(x == cummin(x))

isequal <- function(a, b) {
    attributes(a) <- NULL
    attributes(b) <- NULL
    if (!isTRUE(all.equal(a, b))) {
        print(a, digits = 10)
        print(b, digits = 10)
        FALSE
    } else
        TRUE
}

has_distribution <- function(args)
    !(is.character(args$distribution) && args$distribution[1L] == "none")

check_distribution_arg <- function(distribution,
    values = c("asymptotic", "approximate", "exact", "none")) {
    if (is.character(distribution)) {
        distribution <- match.arg(distribution, values)
        if (distribution == "none")
            function(object) new("NullDistribution")
        else
            eval(call(distribution))
    } else
        distribution
}

setup_args <- function(...) {
    cl <- sys.call(sys.parent())
    fun <- if (inherits(cl$object, "SymmetryProblem"))
               symmetry_test.SymmetryProblem
           else
               independence_test.IndependenceProblem
    cl <- match.call(fun, call = cl, expand.dots = FALSE)
    ## get default arguments and values
    args <- formals(fun)
    args$object <- args$... <- NULL
    nm <- names(args)
    ## replace default values with user-specified values
    for (i in nm[nm %in% names(cl)])
        args[[i]] <- cl[[i]]
    ## override default and user-specified values
    dots <- list(...)
    for (i in nm[nm %in% names(dots)])
        args[[i]] <- dots[[i]]
    lapply(args, eval.parent)
}

statnames <- function(object) {
    xtrans <- object@xtrans
    ytrans <- object@ytrans
    block  <- object@block

    nr <- ncol(xtrans)
    nc <- ncol(ytrans)

    dn <- if (nlevels(block) == 1)
              list(colnames(xtrans), colnames(ytrans))
          else
              list(colnames(xtrans), colnames(ytrans), levels(block))
    if (is.null(dn[[1L]])) {
        if (nr == 1L) {
            dn[[1L]] <- ""
        } else {
            dn[[1L]] <- paste0("X", seq_len(nr))
        }
    }
    if (is.null(dn[[2L]])) {
        if (nc == 1L) {
            dn[[2L]] <- ""
        } else {
            dn[[2L]] <- paste0("Y", seq_len(nc))
        }
    }
    list(dimnames = dn,
         names = paste(rep.int((dn[[1L]]), nc),
                       rep.int((dn[[2L]]), rep.int(nr, nc)),
                       sep = if (all(dn[[1L]] == "") | all(dn[[2L]] == "")) ""
                             else ":"))
}

varnames <- function(object) {
    yordered <- vapply(object@y, is.ordered, NA)
    ynames <- paste0(colnames(object@y), ifelse(yordered, " (ordered)", ""),
                     collapse = ", ")

    if (length(object@x) == 1) {
        if (is.ordered(object@x[[1]])) {
            xnames <- paste0(
                colnames(object@x), " (",
                paste0(levels(object@x[[1]]), collapse = " < "),
                ")"
            )
        } else {
            if (is.factor(object@x[[1]])) {
                xnames <- paste0(
                    colnames(object@x), " (",
                    paste0(levels(object@x[[1]]), collapse = ", "),
                    ")"
                )
            } else {
                xnames <- colnames(object@x)
            }
        }
    } else {
        xordered <- vapply(object@x, is.ordered, NA)
        xnames <- paste0(colnames(object@x), ifelse(xordered, "(ordered)", ""),
                         collapse = ", ")
    }

    if (nlevels(object@block) > 1) {
        bn <- attr(object@block, "blockname")
        if (is.null(bn))
            bn <- "block"
        xnames <- paste(xnames, paste("\n\t stratified by", bn))
    }

    if (nchar(xnames) > getOption("width") / 2)
        paste(ynames, "by\n\t", xnames, collapse = "")
    else
        paste(ynames, "by", xnames, collapse = "")
}

`%EQ%` <- function(x, y)
    abs(x - y) <= sqrt_eps

`%NE%` <- function(x, y)
    abs(x - y) > sqrt_eps

`%GE%` <- function(x, y)
    (y - x) <= sqrt_eps

`%LE%` <- function(x, y)
    (x - y) <= sqrt_eps

`%GT%` <- function(x, y)
    (x - y) > sqrt_eps

`%LT%` <- function(x, y)
    (y - x) > sqrt_eps

### don't use! never!
get_weights <- function(object) object@statistic@weights
get_xtrans <- function(object) object@statistic@xtrans
get_ytrans <- function(object) object@statistic@ytrans

is_unity <- function(x)
    max(abs(x - 1.0)) < sqrt_eps

setRownames <- function(object, value) {
    rownames(object) <- value
    object
}

setColnames <- function(object, value) {
    colnames(object) <- value
    object
}

setDimnames <- function(object, value) {
    dimnames(object) <- value
    object
}

setAttributes <- function(object, value) {
    attributes(object) <- value
    object
}

### heuristic for determining the printed number of decimal digits
### note that, e.g., 1.00 --> 0, 1.10 --> 1, 1.01 --> 2
n_decimal_digits <-
    function(x)
{
    nchar(sub("^-?[[:space:]]?[[:digit:]]*[.]?", "",
              format(x, digits = 15, scientific = FALSE)[1]))
}

### Back-compatibility
if (getRversion() < "4.1.0") {
    ...names <- function() {
        eval(quote(names(list(...))), envir = parent.frame())
    }
}
