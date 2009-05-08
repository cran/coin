

### single step maxT multiple testing procedure
singlestep <- function(object, ...) {

    if (object@statistic@alternative == "two.sided") {
        ts <- abs(statistic(object, "standardized"))
    } else {
        ts <- statistic(object, "standardized")
    }
    ret <- 1 - sapply(ts, pperm, object = object, ...)  
    ret <- matrix(ret, nrow = nrow(ts), ncol = ncol(ts))
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}

### algorithm 2.8 (Free Step-Down Resampling Method) in
### Westfall & Young (1993), page 66 _using standardized 
### statistics instead of p-values_!
### <FIXME>
sdmaxT <- function(pls, ts) {

    ### order of original statistics
    rts <- order(ts)

    ### algorithm 2.8 (Free Step-Down Resampling Method) in 
    ### Westfall & Young (1993), page 66 _using standardized
    ### statistics instead of p-values_!
    q <- pls[,rts, drop = FALSE]

    if (ncol(q) > 1) {
        for (j in 2:ncol(q))
            q[,j] <- pmax(q[,j], q[,j-1])
    }
    ret <- matrix(rowMeans(GE(t(q), ts[rts]))[rank(ts)],
                  nrow = nrow(ts), ncol = ncol(ts))
    
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}
    

### stepdown maxT multiple testing procedure
stepdown <- function(object, ...) {

    if (!(extends(class(object), "MaxTypeIndependenceTest") &&
          extends(class(object@distribution), "ApproxNullDistribution")))
        stop(sQuote("object"), " is not of class ",
             sQuote("MaxTypeIndependenceTest"), 
             " or distribution was not approximated via Monte-Carlo") 
 
    ### raw simulation results, scores have been handled already
    pls <- support(object, raw = TRUE)

    ### standardize
    dcov <- sqrt(variance(object))
    expect <- expectation(object) 
    pls <- t((pls - expect) / dcov)
    ts <- statistic(object, "standardized")

    if (object@statistic@alternative == "two.sided") {
        pls <- abs(pls)
        ts <- abs(ts)
    } 
    if (object@statistic@alternative == "less") {
        pls <- -pls
        ts <- -ts
    }

    sdmaxT(pls, ts)
}

### Bonferroni permutation method (Westfall & Wolfinger, 1997, AmStat 51, 3-8)
dbonf <- function(object, ...) {

   ### <FIXME> this should be possible when the _exact_ marginal
   ### distributions are available
   ### </FIXME>

   if (!(extends(class(object), "MaxTypeIndependenceTest") &&
         extends(class(object@distribution), "ApproxNullDistribution")))
       stop(sQuote("object"), " is not of class ",
            sQuote("MaxTypeIndependenceTest"),
            " or distribution was not approximated via Monte-Carlo")

   alternative <- object@statistic@alternative

   ### standardize
   dcov <- sqrt(variance(object))
   expect <- expectation(object)

   ### raw simulation results, scores have been handled already
   pls <- support(object, raw = TRUE)
   pls <- (pls - expect) / dcov
   ts <- (statistic(object, "standardized"))
   
   pvals <- switch(alternative,
           "less" = rowMeans(LE(pls, drop(ts))),
           "greater" = rowMeans(GE(pls, drop(ts))),
           "two.sided" = rowMeans(GE(abs(pls), abs(drop(ts)))))

   foo <- function(x, t)
       switch(alternative,
           "less" = mean(LE(x, t)),
           "greater" = mean(GE(x, t)),
           "two.sided" = mean(GE(abs(x), abs(t))))

   p <- vector(mode = "list", length = nrow(pls))
   for (i in 1:nrow(pls)) {
       ux <- unique(pls[i,])
       p[[i]] <- sapply(ux, foo, x = pls[i,])
   } 

   ### Bonferroni adjustment (Westfall & Wolfinger, 1997)
   adjp <- rep(1, length(ts))
   for (i in 1:length(pvals)) {
       for (q in 1:length(p)) {
           x <- p[[q]][p[[q]] <= pvals[i]]
           if (length(x) > 0)
               adjp[i] <- adjp[i] * (1 - max(x))
       }
   }
   ret <- matrix(1 - pmin(adjp, 1), nrow = nrow(ts), ncol = ncol(ts))
   rownames(ret) <- rownames(ts)
   colnames(ret) <- colnames(ts)
   ret
}

