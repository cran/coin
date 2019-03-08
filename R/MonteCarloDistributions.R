split_index <-
    function(n, by)
{
    if (n < by)
        by <- n
    lengths(lapply(seq_len(by), function(i) seq.int(i, n, by)),
            use.names = FALSE)
}

MonteCarlo <-
    function(x, y, block, weights, nresample, parallel, ncpus, cl)
{
    montecarlo <- function(nresample)
        .Call(R_PermutedLinearStatistic,
              x, y, weights, integer(0), block, as.double(nresample))

    if (parallel == "no")
        montecarlo(nresample)
    else {
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
            ## advance stream in master process upon exit
            on.exit(.GlobalEnv[[".Random.seed"]] <-
                        nextRNGStream(.GlobalEnv[[".Random.seed"]]))

        if (parallel == "multicore") {
            if (.Platform$OS.type == "windows")
                stop(sQuote(paste0("parallel = ", dQuote("multicore"))),
                     " is not available for MS Windows")
            if (as.integer(ncpus) < 2L)
                warning("parallel operation requires at least two processes")
            do.call("cbind", mclapply(split_index(nresample, ncpus),
                                      FUN = montecarlo, mc.cores = ncpus))
        } else {
            if (is.null(cl)) {
                ## has a default cluster been registered?
                ## see parallel:::defaultCluster
                ## <FIXME> R-3.5.0 has 'getDefaultCluster()' </FIXME>
                cl <- getNamespace("parallel")$.reg$default
                if (is.null(cl)) {
                    ## no default cluster, so setup a PSOCK cluster
                    cl <- makePSOCKcluster(ncpus)
                    on.exit(stopCluster(cl), add = TRUE) # clean-up
                }
            }
            if (RNGkind()[1L] == "L'Ecuyer-CMRG")
                ## distribute streams (using master process) for reproducibility
                clusterSetRNGStream(cl)
            ncpus <- as.integer(length(cl))
            if (ncpus < 2L)
                warning("parallel operation requires at least two processes")
            do.call("cbind", clusterApply(cl, x = split_index(nresample, ncpus),
                                          fun = montecarlo))
        }
    }
}
