
### Conditional Expectation and Covariance
setClass(Class = "ExpectCovar",
    representation = representation(
        expectation = "numeric",
        covariance  = "matrix",
        dimension   = "integer"
   )
)

### Expectation and Covariance of the influence function
### (+ sum of weights)
setClass(Class = "ExpectCovarInfluence",
    representation = representation(
        sumweights = "numeric"
    ),
    contains = "ExpectCovar"
)

### some code from the modeltools package (in preparation) is needed
### modeltools has
# Package: modeltools
# Title: Tools and Classes for Statistical Models
# Date: $Date: 2005/02/21 10:13:12 $
# Version: 0.0-3
# Author: Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> and
#   Friedrich Leisch <Friedrich.Leisch@ci.tuwien.ac.at>
# Maintainer: Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>
# Description: A collection of tools to deal with abstract models
# Depends: R (>= 1.9.0)
# License: GPL

setClass("StatModelCapabilities",
         representation(weights = "logical",
                        subset = "logical"),
         prototype(weights=TRUE, subset=TRUE))

setClass("StatModel",
         representation(name = "character",
                        dpp = "function",
                        fit = "function",
                        capabilities = "StatModelCapabilities"))

## Definition and construction method for class ModelEnv

setClass("FormulaParts",
         representation(formula="list"))
       
setClass("ModelEnv",
         representation(env="environment",
                        get="function",
                        set="function"))

setClass("ModelEnvFormula",
         contains=c("ModelEnv", "FormulaParts"))

setMethod("initialize", "ModelEnv",
function(.Object){

    .Object@env <- new.env()
    
    .Object@get <-
        function(which) get(which, envir=.Object@env, inherits=FALSE)
    
    .Object
})

ModelEnvFormula <- function(formula, data = list(), subset = NULL, 
                            na.action = NULL, frame = NULL,
                            other = list(), designMatrix = TRUE,
                            responseMatrix = TRUE, ...) {
  
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    if (is.null(subset)) mf$subset <- NULL

    ### NA-handling will for the ModelFrame objects later on...
    mf$na.action = na.pass

    MEF = new("ModelEnvFormula")
    MEF@formula <- c(ParseFormula(formula, data=data)@formula, other)
    
    if (is.null(frame)) frame = parent.frame()

    MEF@get <- function(which, data=NULL, frame=parent.frame())
    {
        if(is.null(data))
            RET = get(which, envir=MEF@env, inherits=FALSE)
        else{
            oldData = get(which, envir = MEF@env, inherits=FALSE)
            mf$data = data
            mf$formula = MEF@formula[[which]]
            RET = eval(mf, frame)
            checkData(oldData, RET)
        }
        return(RET)
    }
    
    MEF@set <- function(which=NULL, data=NULL, frame=parent.frame())
    {
        if(is.null(which)) which <- names(MEF@formula)
        if(any(duplicated(which)))
            stop("Some model terms used more than once")
        
        for(name in which){
            
            if(length(MEF@formula[[name]])!=2)
                stop("Invalid formula for ", sQuote(name))
            
            mf$data = data
            mf$formula = MEF@formula[[name]]
            MF = eval(mf, frame)
            if(exists(name, envir=MEF@env, inherits=FALSE))
                checkData(get(name, envir=MEF@env, inherits=FALSE), MF)
            assign(name, MF, envir=MEF@env)
            
            ## <FIXME>
            ## maybe we don't want to save input and response
            ## in the cases below?
            ## </FIXME>
            if(name=="input" && designMatrix){
            assign("designMatrix",
                   model.matrix(attr(MF, "terms"), data=MF),
                   envir=MEF@env)
        }
            if(name=="response" && responseMatrix){
                mt <- attr(MF, "terms")
                attr(mt, "intercept") <- 0
                assign("responseMatrix",
                       model.matrix(mt, data=MF),
                       envir=MEF@env)
            }
        }
    }
    
    MEF@set(which=NULL, data=data, frame=frame)
    
    rm(data)
    rm(other)
    
    ### handle NA's <FIXME>
#    if (!is.null(na.action))
#        MEF = na.action(MEF)
    ## </FIXME>
    MEF
}

checkData <- function(old, new){

  if (!is.null(old)){
    if(!identical(lapply(old, class), lapply(new, class))){
      stop("Classes of new data do not match original data")
    }
    if(!identical(lapply(old, levels), lapply(new, levels))){
      stop("Levels in factors of new data do not match original data")
    }
  }
}
  

ParseFormula <- function(formula, data=list())
{
  formula = terms(formula, data = data)
  attributes(formula) = NULL
  if (length(formula) == 3) {
    fresponse = formula[c(1,2)]
    frhs = formula[c(1,3)]
    if (frhs[[2]] == "1")
      frhs = NULL
  }
  if (length(formula) == 2) {
    fresponse = NULL   
    frhs = formula
  }
  finput = frhs
  fblocks = frhs


  ## <FIXME>
  ## will fail for `y ~ . | blocks' constructs
  ## </FIXME>

  if (!is.null(frhs) && length(frhs[[2]]) > 1) {
    if (deparse(frhs[[2]][[1]]) == "|") {
      finput[[2]] = frhs[[2]][[2]]
      fblocks[[2]] = frhs[[2]][[3]]
    } else {
      fblocks = NULL
    }
  } else {
    fblocks = NULL
  }
  
  fcensored = NULL
  
  if (!is.null(fresponse) && length(fresponse[[2]]) == 3) {
    if (fresponse[[2]][[1]] == "Surv") {
      fcensored = formula(paste("~", fresponse[[2]][[3]]))
      fresponse = formula(paste("~", fresponse[[2]][[2]])) 
    }
  }

  DP = new("FormulaParts")
  
  DP@formula$response=fresponse
  DP@formula$input=finput
  DP@formula$censored=fcensored
  DP@formula$blocks=fblocks

  DP
}

###**********************************************************

setMethod("show", "ModelEnv",
function(object){
    cat("\n")
    cat("A", class(object), "with \n\n")
    n = NULL
    if (has(object, "response")) {
        cat("  response variable(s):  ",
            colnames(object@get("response")), "\n")
        n = nrow(object@get("response"))
    }
    else if (has(object, "responseMatrix")) {
        cat("  response matrix column(s): ",
            colnames(object@get("responseMatrix")), "\n")
        n = nrow(object@get("responseMatrix"))
    }
    
    if (has(object, "input")) {
        cat("  input variable(s):     ",
            colnames(object@get("input")), "\n")
        n = nrow(object@get("input"))
    }
    else if (has(object, "designMatrix")) {
        cat("  design matrix column(s): ",
            colnames(object@get("input")), "\n")
        n = nrow(object@get("designMatrix"))
    }
    
    if (is.null(n)) 
        cat("  no observations\n")
    else
        cat("  number of observations:", n, "\n")
    cat("\n")
})


###**********************************************************

## utility methods for ModelEnv onjects

setGeneric("has", function(object, which) standardGeneric("has"))

setMethod("has", signature(object="ModelEnv", which="character"),
function(object, which){
  exists(which, envir=object@env, inherits=FALSE)
})

setGeneric("dimension", function(object, which) standardGeneric("dimension"))

setMethod("dimension", signature(object="ModelEnv", which="character"),
function(object, which){
  if(has(object, which))
    eval(parse(text=paste("dim(",which,")")) , envir=object@env)
  else
    NULL
})

setGeneric("clone", function(object, ...) standardGeneric("clone"))

## the set() method of ModelEnvFormula objects uses lexical scope on
## various bits and pieces, hence cloning currently returns only a
## ModelEnv object, which only has a trivial get method and no set
## method

setMethod("clone", "ModelEnvFormula",
function(object, copydata=TRUE){

    z <- new("ModelEnv")
    if(copydata){
        for(name in ls(object@env))
            assign(name, object@get(name), envir=z@env)
    }
    z
})

setMethod("subset", "ModelEnvFormula",
function(x, subset, clone=TRUE, ...)
{
    if(clone)
        z <- clone(x, copydata=FALSE)

    for(name in ls(x@env)){
        if(is(x@get(name), "matrix")){
            assign(name, x@get(name)[subset,,drop=FALSE], envir=z@env)
        }
        else{
            assign(name, subset(x@get(name), subset, ...), envir=z@env)
        }
    }
    if(!clone)
        invisible(z)
    else
        return(z)
    
})


###**********************************************************

