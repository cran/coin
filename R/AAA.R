
.onLoad <- function(lib, pkg) {
    .Call("coin_init", PACKAGE = "coin")
    return(TRUE)
}
