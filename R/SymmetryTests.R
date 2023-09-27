### Sign test
sign_test <- function(object, ...) UseMethod("sign_test")

sign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    object <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(object$block)) {
        y <- object$y[[1]]
        if (!is.numeric(y) || is.Surv(y))
            stop(sQuote(colnames(object$y)), " is not a numeric variable")
        x <- object$x[[1]]
        if (!is.numeric(x) || is.Surv(x))
            stop(sQuote(colnames(object$x)), " is not a numeric variable")
        n <- length(x)
        object <- list(x = data.frame(x = gl(2, n)),
                       y = data.frame(y = c(y, x)),
                       block = gl(n, 1, 2 * n))
    }
    object <- new("SymmetryProblem", x = object$x, y = object$y,
                  block = object$block)
    do.call(sign_test, c(object = object, list(...)))
}

sign_test.SymmetryProblem <- function(object, ...) {

    if (!is_numeric_y(object))
        stop(sQuote(colnames(object@y)), " is not a numeric variable")
    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent a paired two-sample problem",
             " (maybe the grouping variable is not a factor?)")

    y <- object@y[[1]]
    x <- object@x[[1]]
    lx <- levels(x); lx1 <- lx[1]; lx2 <- lx[2]
    diffs <- tapply(seq_along(y), object@block, function(b) {
        yb <- y[b]; xb <- x[b]
        yb[xb == lx1] - yb[xb == lx2]
    })
    abs_diffs <- abs(diffs)
    if (all(abs_diffs < eps))
        stop("all pairwise differences equal zero")

    diffs <- diffs[abs_diffs > 0]
    n <- length(diffs)

    object <- new("SymmetryProblem",
                  x = data.frame(x = gl(2, 1, 2 * n, labels = c("pos", "neg"))),
                  y = data.frame(y = as.numeric(rbind(diffs > 0, diffs < 0))),
                  block = gl(n, 2))

    args <- setup_args(teststat = "scalar", paired = TRUE)

    object <- do.call(symmetry_test, c(object = object, args))

    object@method <- "Sign Test"
    object@nullvalue <- 0

    object
}


### Wilcoxon signed-rank test
wilcoxsign_test <- function(object, ...) UseMethod("wilcoxsign_test")

wilcoxsign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    object <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(object$block)) {
        y <- object$y[[1]]
        if (!is.numeric(y) || is.Surv(y))
            stop(sQuote(colnames(object$y)), " is not a numeric variable")
        x <- object$x[[1]]
        if (!is.numeric(x) || is.Surv(x))
            stop(sQuote(colnames(object$x)), " is not a numeric variable")
        n <- length(x)
        object <- list(x = data.frame(x = gl(2, n)),
                       y = data.frame(y = c(y, x)),
                       block = gl(n, 1, 2 * n))
    }
    object <- new("SymmetryProblem", x = object$x, y = object$y,
                  block = object$block)
    do.call(wilcoxsign_test, c(object = object, list(...)))
}

wilcoxsign_test.SymmetryProblem <- function(object,
    zero.method = c("Pratt", "Wilcoxon"), ...) {

    zero.method <- match.arg(zero.method)

    if (!is_numeric_y(object))
        stop(sQuote(colnames(object@y)), " is not a numeric variable")
    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent a paired two-sample problem",
             " (maybe the grouping variable is not a factor?)")

    y <- object@y[[1]]
    x <- object@x[[1]]
    lx <- levels(x); lx1 <- lx[1]; lx2 <- lx[2]
    diffs <- tapply(seq_along(y), object@block, function(b) {
        yb <- y[b]; xb <- x[b]
        yb[xb == lx1] - yb[xb == lx2]
    })
    abs_diffs <- abs(diffs)
    if (all(abs_diffs < eps))
        stop("all pairwise differences equal zero")

    pos_abs_diffs <- abs_diffs > 0
    diffs <- diffs[pos_abs_diffs]
    if (zero.method == "Pratt") {
        rank_abs_diffs <- rank_trafo(abs_diffs)[pos_abs_diffs]
    } else {
        abs_diffs <- abs_diffs[pos_abs_diffs]
        rank_abs_diffs <- rank_trafo(abs_diffs)
    }
    pos <- rank_abs_diffs * (diffs > 0)
    neg <- rank_abs_diffs * (diffs < 0)
    n <- length(pos)

    object <- new("SymmetryProblem",
                  x = data.frame(x = gl(2, 1, 2 * n, labels = c("pos", "neg"))),
                  y = data.frame(y = as.vector(rbind(pos, neg))),
                  block = gl(n, 2))

    args <- setup_args(teststat = "scalar", paired = TRUE)

    object <- do.call(symmetry_test, c(object = object, args))

    if (zero.method == "Pratt")
        object@method <- "Wilcoxon-Pratt Signed-Rank Test"
    else
        object@method <- "Wilcoxon Signed-Rank Test"
    object@nullvalue <- 0

    object
}


### Friedman test
friedman_test <- function(object, ...) UseMethod("friedman_test")

friedman_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("friedman_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

friedman_test.SymmetryProblem <- function(object, ...) {

    block <- object@block
    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = rank_trafo, block = block),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            TRUE
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <- if (is_ordered_x(object)) "scalar"
                     else "quadratic"

    object <- do.call(symmetry_test, c(object = object, args))

    if (is_ordered_x(object@statistic))
        object@method <- "Page Test"
    else
        object@method <- "Friedman Test"

    object
}


### Quade test
quade_test <- function(object, ...) UseMethod("quade_test")

quade_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("quade_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

quade_test.SymmetryProblem <- function(object, ...) {

    block <- object@block
    args <- setup_args(
        ytrafo = function(data) {
            trafo(data, numeric_trafo = function(y) {
                y <- split(y, block)
                R <- lapply(y, function(y) rank_trafo(y) - (length(y) + 1) / 2)
                Q <- rank(vapply(y, function(y) max(y) - min(y), NA_real_,
                                 USE.NAMES = FALSE))
                unsplit(lapply(seq_along(Q), function(i) Q[i] * R[[i]]), block)
            })
        },
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            TRUE
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <- if (is_ordered_x(object)) "scalar"
                     else "quadratic"

    object <- do.call(symmetry_test, c(object = object, args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else
        object@method <- "Quade Test"

    object
}
