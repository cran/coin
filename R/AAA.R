
### attach packages `survival' and `mvtnorm'

.onLoad <- function(lib, pkg) {
    if (!require("survival")) 
        stop("cannot load ", sQuote("survival"))
    if (!require("mvtnorm")) 
        stop("cannot load ", sQuote("mvtnorm"))
    .Call("coin_init", PACKAGE = "coin")
    return(TRUE)
}
