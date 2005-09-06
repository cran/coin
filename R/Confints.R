
confint_location <- function(object, nulldistr, level = 0.95, 
                             approx = FALSE, ...) {

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
        stop("Argument ", sQuote("object"), " is not of class ",
             sQuote("ScalarIndependenceTestStatistic"))

    ### <FIXME> drop unused levels!
    if (!is_2sample(object)) 
        warning(sQuote("object"), " does not represent a two sample problem")
    ### </FIXME>

    if (!extends(class(nulldistr), "NullDistribution"))
        stop("Argument ", sQuote("nulldistr"), " is not of class ",
             sQuote("NullDistribution"))

    if (nlevels(object@block) != 1 || max(abs(object@weights - 1)) > 0)
        stop("cannot compute confidence intervals with blocks or weights")

    alternative <- object@alternative

    if(!((length(level) == 1)         
       && is.finite(level)           
       && (level > 0)          
       && (level < 1)))          
       stop("level must be a single number between 0 and 1")         

    scores <- object@y[[1]]
    groups <- object@xtrans[,1]

    ### raw data
    x <- sort(scores[groups > 0])
    y <- sort(scores[groups < 1])
    alpha <- 1 - level
 
    foo <- function(x, d) x - d

    if (!approx) {

        steps <- c()
        ### explicitely compute all possible steps
        for (lev in levels(object@block)) {
            thisblock <- (object@block == lev)
            ytmp <- sort(split(scores[thisblock], groups[thisblock])[[1]])             
            xtmp <- sort(split(scores[thisblock], groups[thisblock])[[2]])             
            steps <- c(steps, as.vector(outer(xtmp, ytmp, foo)))
        }
        steps <- sort(steps)

        ### computes the statistic under the alternative `d'
        fse <- function(d)
            sum(object@ytrafo(data.frame(c(foo(x,d),y)))[seq(along = x)])

        ### we need to compute the statistics just to the right of
        ### each step       
        ds <- diff(steps)
        justright <- min(abs(ds[ds > .Machine$double.eps]))/2
        jumps <- sapply(steps + justright, fse)

        ### determine if the statistics are in- or decreasing 
        ### jumpsdiffs <- diff(jumps)
        increasing <- all(diff(jumps[c(1,length(jumps))]) > 0)       
        decreasing <- all(diff(jumps[c(1,length(jumps))]) < 0)

        ### this is safe
        if (!(increasing || decreasing)) 
            stop("cannot compute confidence intervals:", 
                 "the step function is not monotone")

        cci <- function(alpha) {
            ### the quantiles:
            ### we reject iff
            ###
            ###   STATISTIC <  qlower OR
            ###   STATISTIC >= qupper
            ###
            qlower <- drop(qperm(nulldistr, alpha/2) * sqrt(variance(object)) +
                           expectation(object))
            qupper <- drop(qperm(nulldistr, 1 - alpha/2) * 
                           sqrt(variance(object)) + expectation(object))

            ## Check if the statistic exceeds both quantiles first.
            if (qlower < min(jumps) || qupper > max(jumps)) {
                warning("Cannot compute confidence intervals")
                return(c(NA, NA))
            }

            if (increasing) {
                ###
                ###  We do NOT reject for all steps with
                ###
                ###     STATISTICS >= qlower AND
                ###     STATISTICS < qupper
                ###
                ###  but the open right interval ends with the
                ###  step with STATISTIC == qupper
                ###
                ci <- c(min(steps[qlower <= jumps]),
                        min(steps[jumps > qupper])) 
            } else {
                ###
                ###  We do NOT reject for all steps with
                ###
                ###     STATISTICS >= qlower AND
                ###     STATISTICS < qupper
                ###
                ###  but the open left interval ends with the
                ###  step with STATISTIC == qupper
                ###
                ci <- c(min(steps[jumps <= qupper]),
                        min(steps[jumps < qlower]))
            }
            ci
        }
        cint <-  switch(alternative, 
                     "two.sided" = cci(alpha),
                     "greater"   = c(cci(alpha*2)[1], Inf), 
                     "less"      = c(-Inf, cci(alpha*2)[2])
                 )
        attr(cint, "conf.level") <- level    

        ### was: median(steps) which will not work for blocks etc.
        u <- jumps - expectation(object)
        sgr <- ifelse(decreasing, min(steps[u <= 0]), max(steps[u <= 0]))
        sle <- ifelse(decreasing, min(steps[u < 0]), min(steps[u > 0]))

        ESTIMATE <- mean(c(sle, sgr), na.rm = TRUE)
        names(ESTIMATE) <- "difference in location"
    } else {
        ### approximate the steps
        ### Here we search the root of the function `fsa' on the set
        ### c(mumin, mumax).
        ##
        ### This returns a value from c(mumin, mumax) for which
        ### the standardized statistic is equal to the
        ### quantile zq.  This means that the statistic is not
        ### within the critical region, and that implies that '
        ### is a confidence limit for the median.

        fsa <- function(d, zq) {
           STAT <- sum(object@ytrafo(data.frame(c(foo(x,d),y)))[seq(along = x)])
           (STAT - expectation(object)) / sqrt(variance(object)) - zq
        }

        mumin <- min(x) - max(y)
        mumax <- max(x) - min(y)
        ccia <- function(alpha) {
            ## Check if the statistic exceeds both quantiles
            ## first: otherwise `uniroot' won't work anyway 
            statu <- fsa(mumin, zq = qperm(nulldistr, alpha/2))
            statl <- fsa(mumax, zq = qperm(nulldistr, 1 - alpha/2))
            if (sign(statu) == sign(statl)) {
                warning(paste("Samples differ in location:",
                              "Cannot compute confidence set,",
                              "returning NA"))
                return(c(NA, NA))
            }
            u <- uniroot(fsa, c(mumin, mumax), 
                         zq = qperm(nulldistr, alpha/2), ...)$root
            l <- uniroot(fsa, c(mumin, mumax), 
                         zq = qperm(nulldistr, 1 - alpha/2), ...)$root
            ## The process of the statistics does not need to be
            ## increasing: sort is ok here.
            sort(c(u, l))
        }
        cint <- switch(alternative, two.sided = {
                    ccia(alpha)
                }, greater= {  
                    c(ccia(alpha*2)[1], Inf)
                }, less= {
                    c(-Inf, ccia(alpha*2)[2])
                })
                attr(cint, "conf.level") <- level
        ## Check if the statistic exceeds both quantiles first.
        statu <- fsa(mumin, zq = 0)
        statl <- fsa(mumax, zq = 0)
        if (sign(statu) == sign(statl)) {
            ESTIMATE <- NA
            warning("Cannot compute estimate, returning NA")
        } else {
            ESTIMATE <- uniroot(fsa, c(mumin, mumax), zq = 0, ...)$root
        }
        names(ESTIMATE) <- "difference in location"
    }

    RET <- list(conf.int = cint, estimate = ESTIMATE)
    RET
}


confint_scale <- function(object, nulldistr, level = 0.95, 
    approx = FALSE, ...) {

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
        stop("Argument ", sQuote("object"), " is not of class ",    
             sQuote("ScalarIndependenceTestStatistic"))     

    if (!extends(class(nulldistr), "NullDistribution"))                               
        stop("Argument ", sQuote("nulldistr"), " is not of class ",
             sQuote("NullDistribution"))                               

    ### <FIXME> drop unused levels!
    if (!is_2sample(object))
        warning(sQuote("object"), " does not represent a two sample problem")
    ### </FIXME>

    if (nlevels(object@block) != 1 || max(abs(object@weights - 1)) > 0)
        stop("cannot compute confidence intervals with blocks or weights")

    alternative <- object@alternative

    if(!((length(level) == 1)         
       && is.finite(level)           
       && (level > 0)          
       && (level < 1)))          
       stop("level must be a single number between 0 and 1")         

    scores <- object@y[[1]]
    groups <- object@xtrans[,1]

    ### raw data
    x <- sort(scores[groups > 0])
    y <- sort(scores[groups < 1])
    alpha <- 1 - level
 
    foo <- function(x, d) x / d

    if (!approx) {

        ### explicitely compute all possible steps
        steps <- c()
        for (lev in levels(object@block)) {
            thisblock <- (object@block == lev)
            ytmp <- sort(split(scores[thisblock], groups[thisblock])[[1]])
            xtmp <- sort(split(scores[thisblock], groups[thisblock])[[2]])
            ratio <-  outer(xtmp, ytmp, "/")
            aratio <- ratio[ratio >= 0]
            steps <- c(steps, aratio)
        }
        steps <- sort(steps)

        ### computes the statistic under the alternative `d'
        fse <- function(d)
            sum(object@ytrafo(data.frame(c(foo(x,d),y)))[seq(along = x)])

        ### we need to compute the statistics just to the right of
        ### each step       
        ds <- diff(steps)
        justright <- min(abs(ds[ds > .Machine$double.eps]))/2
        jumps <- sapply(steps + justright, fse)

        ### determine if the statistics are in- or decreasing 
        ### jumpsdiffs <- diff(jumps)
        increasing <- all(diff(jumps[c(1,length(jumps))]) > 0)       
        decreasing <- all(diff(jumps[c(1,length(jumps))]) < 0)

        ### this is safe
        if (!(increasing || decreasing)) 
            stop("cannot compute confidence intervals:", 
                 "the step function is not monotone")

        cci <- function(alpha) {
            ### the quantiles:
            ### we reject iff
            ###
            ###   STATISTIC <  qlower OR
            ###   STATISTIC >= qupper
            ###
            qlower <- drop(qperm(nulldistr, alpha/2) * sqrt(variance(object)) +
                           expectation(object))
            qupper <- drop(qperm(nulldistr, 1 - alpha/2) * 
                           sqrt(variance(object)) + expectation(object))

            ## Check if the statistic exceeds both quantiles first.
            if (qlower < min(jumps) || qupper > max(jumps)) {
                warning("Cannot compute confidence intervals")
                return(c(NA, NA))
            }

            if (increasing) {
                ###
                ###  We do NOT reject for all steps with
                ###
                ###     STATISTICS >= qlower AND
                ###     STATISTICS < qupper
                ###
                ###  but the open right interval ends with the
                ###  step with STATISTIC == qupper
                ###
                ci <- c(min(steps[qlower <= jumps]), 
                        min(steps[jumps > qupper]))
            } else {
                ###
                ###  We do NOT reject for all steps with
                ###
                ###     STATISTICS >= qlower AND
                ###     STATISTICS < qupper
                ###
                ###  but the open left interval ends with the
                ###  step with STATISTIC == qupper
                ###
                ci <- c(min(steps[jumps <= qupper]),
                        min(steps[jumps < qlower]))
            }
            ci
        }
        cint <-  switch(alternative, 
                     "two.sided" = cci(alpha),
                     "greater"   = c(cci(alpha*2)[1], Inf), 
                     "less"      = c(0, cci(alpha*2)[2])
                 )
        attr(cint, "conf.level") <- level    
        u <- jumps - expectation(object)
        sgr <- ifelse(decreasing, min(steps[u <= 0]), max(steps[u <= 0]))
        sle <- ifelse(decreasing, min(steps[u < 0]), min(steps[u > 0]))

        ESTIMATE <- mean(c(sle, sgr), na.rm = TRUE)
        names(ESTIMATE) <- "ratio of scales"
    } else {
        ### approximate the steps
        ### Here we search the root of the function `fsa' on the set
        ### c(mumin, mumax).
        ##
        ### This returns a value from c(mumin, mumax) for which
        ### the standardized statistic is equal to the
        ### quantile zq.  This means that the statistic is not
        ### within the critical region, and that implies that '
        ### is a confidence limit for the median.

        fsa <- function(d, zq) {
           STAT <- sum(object@ytrafo(data.frame(c(foo(x,d),y)))[seq(along = x)])
           (STAT - expectation(object)) / sqrt(variance(object)) - zq
        }

        srangepos <- NULL
        srangeneg <- NULL
        if (any(x[x > 0]) && any(y[y > 0]))
            srangepos <-
                c(min(x[x>0], na.rm=TRUE)/max(y[y>0], na.rm=TRUE),
                  max(x[x>0], na.rm=TRUE)/min(y[y>0], na.rm=TRUE))
        if (any(x[x <= 0]) && any(y[y < 0]))
            srangeneg <-
                c(min(x[x<=0], na.rm=TRUE)/max(y[y<0], na.rm=TRUE),
                  max(x[x<=0], na.rm=TRUE)/min(y[y<0], na.rm=TRUE))
        if (any(is.infinite(c(srangepos, srangeneg)))) {
            stop(paste("Cannot compute asymptotic confidence",
                       "set or estimator"))
        }
        mumin <- range(c(srangepos, srangeneg), na.rm=FALSE)[1]
        mumax <- range(c(srangepos, srangeneg), na.rm=FALSE)[2]

        ccia <- function(alpha) {
            ## Check if the statistic exceeds both quantiles
            ## first: otherwise `uniroot' won't work anyway 
            statu <- fsa(mumin, zq = qperm(nulldistr, alpha/2))
            statl <- fsa(mumax, zq = qperm(nulldistr, 1 - alpha/2))
            if (sign(statu) == sign(statl)) {
                warning(paste("Samples differ in location:",
                              "Cannot compute confidence set,",
                              "returning NA"))
                return(c(NA, NA))
            }
            u <- uniroot(fsa, c(mumin, mumax), 
                         zq = qperm(nulldistr, alpha/2), ...)$root
            l <- uniroot(fsa, c(mumin, mumax), 
                         zq = qperm(nulldistr, 1 - alpha/2), ...)$root
            ## The process of the statistics does not need to be
            ## increasing: sort is ok here.
            sort(c(u, l))
        }
        cint <- switch(alternative, two.sided = {
                    ccia(alpha)
                }, greater= {  
                    c(ccia(alpha*2)[1], Inf)
                }, less= {
                    c(0, ccia(alpha*2)[2])
                })
                attr(cint, "conf.level") <- level
        ## Check if the statistic exceeds both quantiles first.
        statu <- fsa(mumin, zq = 0)
        statl <- fsa(mumax, zq = 0)
        if (sign(statu) == sign(statl)) {
            ESTIMATE <- NA
            warning("Cannot compute estimate, returning NA")
        } else {
            ESTIMATE <- uniroot(fsa, c(mumin, mumax), zq = 0, ...)$root
        }
        names(ESTIMATE) <- "ratio of scales"
    }

    RET <- list(conf.int = cint, estimate = ESTIMATE)
    RET
}

simconfint_location <- function(object, level = 0.95, 
    approx = FALSE, ...) {

    if (!(is_Ksample(object@statistic) && 
        extends(class(object), "MaxTypeIndependenceTest")))
        stop(sQuote("object"), " is not an object of class ",
             sQuote("MaxTypeIndependenceTest"), 
             " representing a K sample problem")

    xtrans <- object@statistic@xtrans
    if (!all(apply(xtrans, 2, function(x) all(x %in% c(-1, 0, 1)))))
        stop("Only differences are allowed as contrasts")

    estimate <- c()
    lower <- c()
    upper <- c()

    ### transform max(abs(x))-type distribution into a 
    ### distribution symmetric around zero
    nnd <- object@distribution
    nnd@q <- function(p) {
        pp <- p
        if (p > 0.5)
            pp <- 1 - pp
        q <- qperm(object@distribution, 1 - pp)
        if (p < 0.5)
        return(-q)
        else
        return(q)
    }
        
    for (i in 1:ncol(xtrans)) {
        thisset <- abs(xtrans[,i]) > 0
        ip <- new("IndependenceProblem", 
                  object@statistic@x[thisset,,drop = FALSE],
                  object@statistic@y[thisset,,drop = FALSE],
                  object@statistic@block[thisset])
        
        itp <- independence_test(ip, teststat = "scalar", 
            distribution = "asympt", alternative = "two.sided", 
            yfun = object@statistic@ytrafo, ...)

        ci <- confint_location(itp@statistic, nnd,
                               level = level, approx =approx, ...)
        estimate <- c(estimate, ci$estimate)
        lower <- c(lower, ci$conf.int[1])
        upper <- c(upper, ci$conf.int[2])
    }
    RET <- data.frame(Estimate = estimate, lower = lower, upper = upper)
    colnames(RET)[2:3] <- paste(c((1 - level)/2, 1 - (1 - level)/2)*100, "%")
    rownames(RET) <- colnames(object@statistic@xtrans)
    attr(RET, "conf.level") <- level
    return(RET)
}

