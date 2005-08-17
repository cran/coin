

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
    if (object@statistic@alternative == "two.sided") {
        pls <- lapply(pls, function(x)
            (abs(x - expect) / dcov)
        )
        ts <- abs(statistic(object, "standardized"))
    } else {
        pls <- lapply(pls, function(x)
            (x - expect) / dcov
        )
        ts <- statistic(object, "standardized")
    }

    ### order of original statistics
    rts <- order(ts)
    if (object@statistic@alternative == "less") rts <- rev(rts)

    ### algorithm 2.8 (Free Step-Down Resampling Method) in
    ### Westfall & Young (1993), page 66 _using standardized 
    ### statistics instead of p-values_!
    q <- matrix(unlist(pls), nrow = length(pls), 
                byrow = TRUE)[,rts]

    if (object@statistic@alternative == "less") {
        for (j in 2:ncol(q))
            q[,j] <- pmin(q[,j], q[,j-1])
        ret <- matrix(rowMeans(t(q) < ts[rts])[rev(rank(ts))], 
                      nrow = nrow(ts), ncol = ncol(ts))
    } else {
        for (j in 2:ncol(q))
            q[,j] <- pmax(q[,j], q[,j-1])
        ret <- matrix(rowMeans(t(q) >= ts[rts])[rank(ts)], 
                      nrow = nrow(ts), ncol = ncol(ts))
    }
    
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
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
   pls <- lapply(pls, function(x)
            (x - expect) / dcov)
   pls <- t(matrix(unlist(pls), nrow = length(pls), byrow = TRUE))
   ts <- (statistic(object, "standardized"))
   
   pvals <- switch(alternative,
           "less" = rowMeans(pls <= drop(ts)),
           "greater" = rowMeans(pls >= drop(ts)),
           "two.sided" = rowMeans(abs(pls) >= abs(drop(ts))))

   foo <- function(x, t)
       switch(alternative,
           "less" = mean(x <= t),
           "greater" = mean(x >= t),
           "two.sided" = mean(abs(x) >= abs(t)))

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

