
### generic method for asymptotic null distributions
setGeneric("AsymptNullDistribution", function(object, ...)
    standardGeneric("AsymptNullDistribution"))

### method for max-type test statistics
setMethod(f = "AsymptNullDistribution", 
          signature = "ScalarIndependenceTestStatistic", 
          definition = function(object, ...) {

              RET <- new("AsymptNullDistribution")
              RET@p <- function(q) pnorm(q)
              RET@q <- function(p) qnorm(p)
              RET@d <- function(x) dnorm(x)
              RET@pvalue <- function(q) {
                  switch(object@alternative,
                      "less"      = pnorm(q),
                      "greater"   = 1 - pnorm(q),
                      "two.sided" = 2 * min(pnorm(q), 1 - pnorm(q))  
                  )
              }
              RET@support <- function(p = 1e-5) c(RET@q(p), RET@q(1 - p))
              RET@name <- "normal distribution"
              return(RET)
          }
)

### just a wrapper
pmv <- function(lower, upper, mean, corr, ...) {
    if (length(corr) > 1) {
        pmvnorm(lower = lower, upper = upper, mean = mean, corr = corr, ...)
    } else {
        pmvnorm(lower = lower, upper = upper, mean = mean, sigma = 1, ...)
    }
} 


### method for max-type test statistics
setMethod(f = "AsymptNullDistribution", 
          signature = "MaxTypeIndependenceTestStatistic", 
          definition = function(object, ...) {

              corr <- cov2cor(covariance(object))
              pq <- length(expectation(object))
              RET <- new("AsymptNullDistribution")
              RET@p <- function(q) {
                  p <- switch(object@alternative,
                      "less"      = pmv(lower = q, upper = Inf, 
                                        mean = rep(0, pq),
                                        corr = corr, ...),
                      "greater"   = pmv(lower = -Inf, upper = q, 
                                        mean = rep(0, pq),
                                        corr = corr, ...),
                      "two.sided" = pmv(lower = -abs(q), upper = abs(q),
                                        mean = rep(0, pq),
                                        corr = corr, ...)
                  )
                  error <- attr(p, "error")
                  attr(p, "error") <- NULL
                  ci <- c(max(0, p - error), min(p + error, 1))
                  attr(ci, "conf.level") <- 0.99
                  attr(p, "conf.int") <- ci
                  class(p) <- "MCp"
                  p
              }
              RET@q <- function(p) {
                  if (length(corr) > 1) 
                      q <- qmvnorm(p, mean = rep(0, pq), 
                              corr = corr, tail = "both.tails", ...)$quantile
                  else
                      q <- qmvnorm(p, mean = rep(0, pq),
                              sigma = 1, tail = "both.tails", ...)$quantile   

                  attributes(q) <- NULL
                  q
              }
              RET@d <- function(x) dmvnorm(x)
              RET@pvalue <- function(q) {
                  p <- 1 - RET@p(q)
                  attr(p, "conf.int") <- 1 - attr(p, "conf.int")[c(2,1)]
                  attr(attr(p, "conf.int"), "conf.level") <- 0.99
                  class(p) <- "MCp"
                  p
              }

              RET@support <- function(p = 1e-5) c(RET@q(p), RET@q(1 - p))

              RET@name <- "multivariate normal distribution"
              RET@parameters <- list(corr = corr)
              return(RET)
          }
)

### method for quad-type test statistics
setMethod(f = "AsymptNullDistribution", 
          signature = "QuadTypeIndependenceTestStatistic", 
          definition = function(object, ...) {

              RET <- new("AsymptNullDistribution")
              RET@p <- function(q) pchisq(q, df = object@df)
              RET@q <- function(p) qchisq(p, df = object@df)
              RET@d <- function(d) dchisq(d, df = object@df)
              RET@pvalue <- function(q) 1 - RET@p(q)

              RET@support <- function(p = 1e-5) c(0, RET@q(1 - p))

              RET@name <- "chisq distribution"
              RET@parameters <- list(df = object@df)
              return(RET)
          }
)

### generic method for exact null distributions
setGeneric("ExactNullDistribution", function(object, ...)
    standardGeneric("ExactNullDistribution"))

setMethod(f = "ExactNullDistribution",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, algorithm = c("shift", "split-up"), 
                                ...) {

              if (is_2sample(object)) {
                  if (algorithm == "shift")
                      return(SR_shift_2sample(object, ...))
                  if (algorithm == "split-up")
                      return(vdW_split_up_2sample(object))
              }
              stop(sQuote("object"), " is not a two sample problem")

          }
)

### generic method for approximate null distributions
setGeneric("ApproxNullDistribution", function(object, ...)
    standardGeneric("ApproxNullDistribution"))

setMethod(f = "ApproxNullDistribution",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, B = 1000, ...) {

              if (!(max(abs(object@weights - 1.0)) < eps()))
                  stop("cannot approximate distribution with non-unity weights")

              pls <- plsraw <- .Call("R_MonteCarloIndependenceTest", object@xtrans, 
                  object@ytrans, as.integer(object@block), as.integer(B), 
                  PACKAGE = "coin")

              ### <FIXME> can transform p, q, x instead of those </FIXME>
              pls <- sort(round((pls - expectation(object)) / 
                         sqrt(variance(object)), 10))

              RET <- new("ApproxNullDistribution")

              RET@p <- function(q) {
                  p <- mean(pls <= round(q, 10))
                  attr(p, "conf.int") <- binom.test(round(p * B), B, 
                      conf.level = 0.99)$conf.int
                  class(p) <- "MCp"
                  p
              }

              RET@q <- function(p) pls[length(pls) * p]
              RET@d <- function(x) {
                  tmp <- abs(pls - x)
                  mean(tmp == tmp[which.min(tmp)])
              }
              RET@pvalue <- function(q) {
                  p <- switch(object@alternative,
                      "less"      = mean(pls <= round(q, 10)), 
                      "greater"   = mean(pls >= round(q, 10)),
                      "two.sided" = mean(abs(pls) >= round(abs(q), 10)))
                  attr(p, "conf.int") <- binom.test(round(p * B), B, 
                      conf.level = 0.99)$conf.int
                  class(p) <- "MCp"
                  p
              }
              RET@support <- function(raw = FALSE) {
                  if (raw) return(plsraw)
                  sort(unique(drop(pls)))
              }
              return(RET)
          }
)

setMethod(f = "ApproxNullDistribution",
          signature = "MaxTypeIndependenceTestStatistic",
          definition = function(object, B = 1000, ...) {

              if (!(max(abs(object@weights - 1.0)) < eps()))
                  stop("cannot approximate distribution with non-unity weights")

              pls <- plsraw <- .Call("R_MonteCarloIndependenceTest", object@xtrans, 
                  object@ytrans, as.integer(object@block), as.integer(B), 
                  PACKAGE = "coin")

              fun <- switch(object@alternative,
                  "less" = min,
                  "greater" = max,
                  "two.sided" = function(x) max(abs(x))
              )

              dcov <- sqrt(variance(object))
              expect <- expectation(object)
              pls <- (pls - expect) / dcov
              pls <- switch(object@alternative,
                  "less" = do.call("pmin", as.data.frame(t(pls))),
                  "greater" = do.call("pmax", as.data.frame(t(pls))),
                  "two.sided" = do.call("pmax", as.data.frame(t(abs(pls)))))
              pls <- sort(round(pls, 10))

              RET <- new("ApproxNullDistribution")

              RET@p <- function(q) {
                  p <- switch(object@alternative,
                      "less" = mean(pls >= round(q, 10)),
                      "greater" = mean(pls <= round(q, 10)),
                      "two.sided" = mean(pls <= round(q, 10))
                  )
                  attr(p, "conf.int") <- binom.test(round(p * B), B, 
                      conf.level = 0.99)$conf.int
                  class(p) <- "MCp"
                  p
              }

              RET@q <- function(p) pls[length(pls) * p]
              RET@d <- function(x) {
                  tmp <- abs(pls - x)
                  mean(tmp == tmp[which.min(tmp)])
              }
              RET@pvalue <- function(q) {
                  p <- switch(object@alternative,
                      "less" = mean(pls <= round(q, 10)),
                      "greater" = mean(pls >= round(q, 10)),
                      "two.sided" = mean(pls >= round(q, 10))
                  )
                  attr(p, "conf.int") <- binom.test(round(p * B), B,
                      conf.level = 0.99)$conf.int
                  class(p) <- "MCp"
                  p
              }

              RET@support <- function(raw = FALSE) {
                  if (raw) return(plsraw)
                  sort(unique(drop(pls)))
              }
              RET@name = "MonteCarlo distribution"
              return(RET)
          }
)

setMethod(f = "ApproxNullDistribution",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, B = 1000, ...) {

              if (!(max(abs(object@weights - 1.0)) < eps()))
                  stop("cannot approximate distribution with non-unity weights")

              pls <- plsraw <- .Call("R_MonteCarloIndependenceTest", object@xtrans, 
                  object@ytrans, as.integer(object@block), as.integer(B), 
                  PACKAGE = "coin")

              dcov <- object@covarianceplus
              expect <- expectation(object)
              a <- pls - expect
              pls <- rowSums((t(a) %*% dcov) * t(a))
              pls <- sort(round(pls, 10))

              RET <- new("ApproxNullDistribution")

              RET@p <- function(q) {
                  p <- mean(pls <= round(q, 10))
                  attr(p, "conf.int") <- binom.test(round(p * B), B, 
                      conf.level = 0.99)$conf.int
                  class(p) <- "MCp"
                  p
              }

              RET@q <- function(p) pls[length(pls) * p]
              RET@d <- function(x) {
                  tmp <- abs(pls - x)
                  mean(tmp == tmp[which.min(tmp)])
              }
              RET@pvalue <- function(q) {
                  p <- mean(pls >= round(q, 10))
                  attr(p, "conf.int") <- binom.test(round(p * B), B, 
                      conf.level = 0.99)$conf.int
                  class(p) <- "MCp"
                  p
              }
              RET@support <- function(raw = FALSE) {
                  if (raw) return(plsraw)
                  sort(unique(drop(pls)))
              }
              RET@name = "MonteCarlo distribution"
              return(RET)
          }
)

confint.ScalarIndependenceTestConfint <- function(object, parm, level = 0.95, 
    ...) {
        if ("level" %in% names(match.call()))
            x <- object@confint(level)
        else
            x <- object@confint(object@conf.level)
        class(x) <- "ci"
        return(x)
}
