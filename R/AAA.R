.onLoad <- function(lib, pkg)
    .Call(coin_init)

.onUnload <- function(libpath)
    library.dynam.unload("coin", libpath)
