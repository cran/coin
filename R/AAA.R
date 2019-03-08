.onLoad <- function(libname, pkgname) {
    ns <- getNamespace(pkgname)
    assign(     "eps",      .Machine$double.eps , envir = ns)
    assign("sqrt_eps", sqrt(.Machine$double.eps), envir = ns)
}

.onUnload <- function(libpath)
    library.dynam.unload("coin", libpath)
