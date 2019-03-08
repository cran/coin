### Regression tests for naming of results

set.seed(290875)
library("coin")

y1 <- sample(1:20)
y2 <- rnorm(20)
x <- gl(2, 10)


###
### Asymptotic
###

### Scalar
asy_scl <- independence_test(y1 ~ x, distr = "asymptotic", teststat = "scalar")
s <- statistic(asy_scl)

stopifnot(identical(   names(    expectation(asy_scl                         )), "1"            ))
stopifnot(identical(   names(       variance(asy_scl                         )), "1"            ))
stopifnot(identical(dimnames(     covariance(asy_scl                         )), list("1", "1") ))
stopifnot(  is.null(   names(      statistic(asy_scl,   type = "test"        ))                 )) # named in < 1.3-0
stopifnot(identical(dimnames(      statistic(asy_scl,   type = "linear"      )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(asy_scl,   type = "centered"    )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(asy_scl,   type = "standardized")), list("1", "")  ))
## stopifnot(  is.null(   names(        support(asy_scl                         ))                 )) # NA
stopifnot(  is.null(   names(          dperm(asy_scl,      x =       s       ))                 ))
stopifnot(identical(   names(          dperm(asy_scl,      x = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          pperm(asy_scl,      q =       s       ))                 ))
stopifnot(identical(   names(          pperm(asy_scl,      q = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          qperm(asy_scl,      p =       0.75    ))                 ))
stopifnot(identical(   names(          qperm(asy_scl,      p = c(s = 0.75)   )), "s"            ))
stopifnot(  is.null(   names(          rperm(asy_scl,      n =       5       ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_scl, method = "global"      ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_scl, method = "single-step" ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_scl, method = "step-down"   ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_scl, method = "unadjusted"  ))                 ))
## stopifnot(  is.null(   names(      midpvalue(asy_scl                         ))                 )) # NA
## stopifnot(identical(   names(pvalue_interval(asy_scl                         )), c("p_0", "p_1"))) # NA
## stopifnot(  is.null(   names(           size(asy_scl,  alpha = 0.05          ))                 )) # NA

### Quadratic univariate
asy_qdr_u <- independence_test(y1 ~ x, distr = "asymptotic", teststat = "quadratic")
s <- statistic(asy_qdr_u)

stopifnot(identical(   names(    expectation(asy_qdr_u                         )), "1"            ))
stopifnot(identical(   names(       variance(asy_qdr_u                         )), "1"            ))
stopifnot(identical(dimnames(     covariance(asy_qdr_u                         )), list("1" , "1")))
stopifnot(  is.null(   names(      statistic(asy_qdr_u,   type = "test"        ))                 ))
stopifnot(identical(dimnames(      statistic(asy_qdr_u,   type = "linear"      )), list("1" , "") ))
stopifnot(identical(dimnames(      statistic(asy_qdr_u,   type = "centered"    )), list("1" , "") ))
stopifnot(identical(dimnames(      statistic(asy_qdr_u,   type = "standardized")), list("1" , "") ))
## stopifnot(  is.null(   names(        support(asy_qdr_u                         ))                 )) # NA
stopifnot(  is.null(   names(          dperm(asy_qdr_u,      x =       s       ))                 ))
stopifnot(identical(   names(          dperm(asy_qdr_u,      x = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          pperm(asy_qdr_u,      q =       s       ))                 ))
stopifnot(identical(   names(          pperm(asy_qdr_u,      q = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          qperm(asy_qdr_u,      p =       0.75    ))                 ))
stopifnot(identical(   names(          qperm(asy_qdr_u,      p = c(s = 0.75)   )), "s"            ))
stopifnot(  is.null(   names(          rperm(asy_qdr_u,      n =       5       ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_u, method = "global"      ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_u, method = "single-step" ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_u, method = "step-down"   ))                 ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_u, method = "unadjusted"  ))                 ))
## stopifnot(  is.null(   names(      midpvalue(asy_qdr_u                         ))                 )) # NA
## stopifnot(identical(   names(pvalue_interval(asy_qdr_u                         )), c("p_0", "p_1"))) # NA
## stopifnot(  is.null(   names(           size(asy_qdr_u,  alpha = 0.05          ))                 )) # NA

### Maximum
asy_mxm <- independence_test(y1 + y2 ~ x, distribution = "asymptotic", teststat = "maximum")
s <- statistic(asy_mxm)

stopifnot(identical(   names(    expectation(asy_mxm                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(   names(       variance(asy_mxm                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(dimnames(     covariance(asy_mxm                         )), list(c("1:y1", "1:y2"), c("1:y1", "1:y2"))))
stopifnot(  is.null(   names(      statistic(asy_mxm,   type = "test"        ))                                            ))
stopifnot(identical(dimnames(      statistic(asy_mxm,   type = "linear"      )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(asy_mxm,   type = "centered"    )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(asy_mxm,   type = "standardized")), list("1", c("y1", "y2"))                  ))
## stopifnot(  is.null(   names(        support(asy_mxm                         ))                                            )) # NA
## stopifnot(  is.null(   names(          dperm(asy_mxm,      x =       s       ))                                            )) # NA
## stopifnot(identical(   names(          dperm(asy_mxm,      x = c(s = s)      )), "s"                                       )) # NA
stopifnot(  is.null(   names(          pperm(asy_mxm,      q =       s       ))                                            ))
stopifnot(identical(   names(          pperm(asy_mxm,      q = c(s = s)      )), "s"                                       ))
stopifnot(  is.null(   names(          qperm(asy_mxm,      p =       0.75    ))                                            ))
stopifnot(identical(   names(          qperm(asy_mxm,      p = c(s = 0.75)   )), "s"                                       ))
## stopifnot(  is.null(   names(          rperm(asy_mxm,      n =       5       ))                                            )) # NA
stopifnot(  is.null(   names(         pvalue(asy_mxm, method = "global"      ))                                            ))
stopifnot(identical(dimnames(         pvalue(asy_mxm, method = "single-step" )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(         pvalue(asy_mxm, method = "step-down"   )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(         pvalue(asy_mxm, method = "unadjusted"  )), list("1", c("y1", "y2"))                  ))
## stopifnot(  is.null(   names(      midpvalue(asy_mxm                         ))                                            )) # NA
## stopifnot(identical(   names(pvalue_interval(asy_mxm                         )), c("p_0", "p_1")                           )) # NA
## stopifnot(  is.null(   names(           size(asy_mxm,  alpha = 0.05          ))                                            )) # NA

### Quadratic multivariate
asy_qdr_m <- independence_test(y1 + y2 ~ x, distr = "asymptotic", teststat = "quadratic")
s <- statistic(asy_qdr_m)

stopifnot(identical(   names(    expectation(asy_qdr_m                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(   names(       variance(asy_qdr_m                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(dimnames(     covariance(asy_qdr_m                         )), list(c("1:y1", "1:y2"), c("1:y1", "1:y2"))))
stopifnot(  is.null(   names(      statistic(asy_qdr_m,   type = "test"        ))                                            ))
stopifnot(identical(dimnames(      statistic(asy_qdr_m,   type = "linear"      )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(asy_qdr_m,   type = "centered"    )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(asy_qdr_m,   type = "standardized")), list("1", c("y1", "y2"))                  ))
## stopifnot(  is.null(   names(        support(asy_qdr_m                         ))                                            )) # NA
stopifnot(  is.null(   names(          dperm(asy_qdr_m,      x =       s       ))                                            ))
stopifnot(identical(   names(          dperm(asy_qdr_m,      x = c(s = s)      )), "s"                                       ))
stopifnot(  is.null(   names(          pperm(asy_qdr_m,      q =       s       ))                                            ))
stopifnot(identical(   names(          pperm(asy_qdr_m,      q = c(s = s)      )), "s"                                       ))
stopifnot(  is.null(   names(          qperm(asy_qdr_m,      p =       0.75    ))                                            ))
stopifnot(identical(   names(          qperm(asy_qdr_m,      p = c(s = 0.75)   )), "s"                                       ))
stopifnot(  is.null(   names(          rperm(asy_qdr_m,      n =       5       ))                                            ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_m, method = "global"      ))                                            ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_m, method = "single-step" ))                                            ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_m, method = "step-down"   ))                                            ))
stopifnot(  is.null(   names(         pvalue(asy_qdr_m, method = "unadjusted"  ))                                            ))
## stopifnot(  is.null(   names(      midpvalue(asy_qdr_m                         ))                                            )) # NA
## stopifnot(identical(   names(pvalue_interval(asy_qdr_m                         )), c("p_0", "p_1" )                          )) # NA
## stopifnot(  is.null(   names(           size(asy_qdr_m,  alpha = 0.05          ))                                            )) # NA


###
### Approximate
###

### Scalar
app_scl <- independence_test(y1 ~ x, distr = "approximate")
s <- statistic(app_scl)

stopifnot(identical(   names(    expectation(app_scl                         )), "1"            ))
stopifnot(identical(   names(       variance(app_scl                         )), "1"            ))
stopifnot(identical(dimnames(     covariance(app_scl                         )), list("1", "1") ))
stopifnot(  is.null(   names(      statistic(app_scl,   type = "test"        ))                 )) # named in < 1.3-0
stopifnot(identical(dimnames(      statistic(app_scl,   type = "linear"      )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(app_scl,   type = "centered"    )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(app_scl,   type = "standardized")), list("1", "")  ))
stopifnot(  is.null(   names(        support(app_scl                         ))                 ))
stopifnot(  is.null(   names(          dperm(app_scl,      x =       s       ))                 ))
stopifnot(identical(   names(          dperm(app_scl,      x = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          pperm(app_scl,      q =       s       ))                 ))
stopifnot(identical(   names(          pperm(app_scl,      q = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          qperm(app_scl,      p =       0.75    ))                 ))
stopifnot(identical(   names(          qperm(app_scl,      p = c(s = 0.75)   )), "s"            )) # unnamed in < 1.3-0
stopifnot(  is.null(   names(          rperm(app_scl,      n =       5       ))                 ))
stopifnot(  is.null(   names(         pvalue(app_scl, method = "global"      ))                 ))
stopifnot(  is.null(   names(         pvalue(app_scl, method = "single-step" ))                 ))
stopifnot(  is.null(   names(         pvalue(app_scl, method = "step-down"   ))                 ))
stopifnot(  is.null(   names(         pvalue(app_scl, method = "unadjusted"  ))                 ))
stopifnot(  is.null(   names(      midpvalue(app_scl                         ))                 ))
stopifnot(identical(   names(pvalue_interval(app_scl                         )), c("p_0", "p_1")))
stopifnot(  is.null(   names(           size(app_scl,  alpha = 0.05          ))                 ))

### Quadratic univariate
app_qdr_u <- independence_test(y1 ~ x, distr = "approximate", teststat = "quadratic")
s <- statistic(app_qdr_u)

stopifnot(identical(   names(    expectation(app_qdr_u                         )), "1"            ))
stopifnot(identical(   names(       variance(app_qdr_u                         )), "1"            ))
stopifnot(identical(dimnames(     covariance(app_qdr_u                         )), list("1", "1") ))
stopifnot(  is.null(   names(      statistic(app_qdr_u,   type = "test"        ))                 ))
stopifnot(identical(dimnames(      statistic(app_qdr_u,   type = "linear"      )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(app_qdr_u,   type = "centered"    )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(app_qdr_u,   type = "standardized")), list("1", "")  ))
stopifnot(  is.null(   names(        support(app_qdr_u                         ))                 ))
stopifnot(  is.null(   names(          dperm(app_qdr_u,      x =       s       ))                 ))
stopifnot(identical(   names(          dperm(app_qdr_u,      x = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          pperm(app_qdr_u,      q =       s       ))                 ))
stopifnot(identical(   names(          pperm(app_qdr_u,      q = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          qperm(app_qdr_u,      p =       0.75    ))                 ))
stopifnot(identical(   names(          qperm(app_qdr_u,      p = c(s = 0.75)   )), "s"            )) # unnamed in < 1.3-0
stopifnot(  is.null(   names(          rperm(app_qdr_u,      n =       5       ))                 ))
stopifnot(  is.null(   names(         pvalue(app_qdr_u, method = "global"      ))                 ))
stopifnot(  is.null(   names(         pvalue(app_qdr_u, method = "single-step" ))                 ))
stopifnot(  is.null(   names(         pvalue(app_qdr_u, method = "step-down"   ))                 ))
stopifnot(  is.null(   names(         pvalue(app_qdr_u, method = "unadjusted"  ))                 ))
stopifnot(  is.null(   names(      midpvalue(app_qdr_u                         ))                 ))
stopifnot(identical(   names(pvalue_interval(app_qdr_u                         )), c("p_0", "p_1")))
stopifnot(  is.null(   names(           size(app_qdr_u,  alpha = 0.05          ))                 ))

### Maximum
app_mxm <- independence_test(y1 + y2 ~ x, distr = "approximate", teststat = "maximum")
s <- statistic(app_mxm)

stopifnot(identical(   names(    expectation(app_mxm                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(   names(       variance(app_mxm                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(dimnames(     covariance(app_mxm                         )), list(c("1:y1", "1:y2"), c("1:y1", "1:y2"))))
stopifnot(  is.null(   names(      statistic(app_mxm,   type = "test"        ))                                            ))
stopifnot(identical(dimnames(      statistic(app_mxm,   type = "linear"      )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(app_mxm,   type = "centered"    )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(app_mxm,   type = "standardized")), list("1", c("y1", "y2"))                  ))
stopifnot(  is.null(   names(        support(app_mxm                         ))                                            ))
stopifnot(  is.null(   names(          dperm(app_mxm,      x =       s       ))                                            ))
stopifnot(identical(   names(          dperm(app_mxm,      x = c(s = s)      )), "s"                                       ))
stopifnot(  is.null(   names(          pperm(app_mxm,      q =       s       ))                                            ))
stopifnot(identical(   names(          pperm(app_mxm,      q = c(s = s)      )), "s"                                       ))
stopifnot(  is.null(   names(          qperm(app_mxm,      p =       0.75    ))                                            ))
stopifnot(identical(   names(          qperm(app_mxm,      p = c(s = 0.75)   )), "s"                                       )) # unnamed in < 1.3-0
stopifnot(  is.null(   names(          rperm(app_mxm,      n =       5       ))                                            ))
stopifnot(  is.null(   names(         pvalue(app_mxm, method = "global"      ))                                            ))
stopifnot(identical(dimnames(         pvalue(app_mxm, method = "single-step" )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(         pvalue(app_mxm, method = "step-down"   )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(         pvalue(app_mxm, method = "unadjusted"  )), list("1", c("y1", "y2"))                  ))
stopifnot(  is.null(   names(      midpvalue(app_mxm                         ))                                            ))
stopifnot(identical(   names(pvalue_interval(app_mxm                         )), c("p_0", "p_1")                           ))
stopifnot(  is.null(   names(           size(app_mxm,  alpha = 0.05          ))                                            ))

### Quadratic multivariate
app_qdr_m <- independence_test(y1 + y2 ~ x, distr = "approximate", teststat = "quadratic")
s <- statistic(app_qdr_m)

stopifnot(identical(   names(    expectation(app_qdr_m                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(   names(       variance(app_qdr_m                         )), c("1:y1", "1:y2")                         ))
stopifnot(identical(dimnames(     covariance(app_qdr_m                         )), list(c("1:y1", "1:y2"), c("1:y1", "1:y2"))))
stopifnot(  is.null(   names(      statistic(app_qdr_m,   type = "test"        ))                                            ))
stopifnot(identical(dimnames(      statistic(app_qdr_m,   type = "linear"      )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(app_qdr_m,   type = "centered"    )), list("1", c("y1", "y2"))                  ))
stopifnot(identical(dimnames(      statistic(app_qdr_m,   type = "standardized")), list("1", c("y1", "y2"))                  ))
stopifnot(  is.null(   names(        support(app_qdr_m                         ))                                            ))
stopifnot(  is.null(   names(          dperm(app_qdr_m,      x =       s       ))                                            ))
stopifnot(identical(   names(          dperm(app_qdr_m,      x = c(s = s)      )), "s"                                       ))
stopifnot(  is.null(   names(          pperm(app_qdr_m,      q =       s       ))                                            ))
stopifnot(identical(   names(          pperm(app_qdr_m,      q = c(s = s)      )), "s"                                       ))
stopifnot(  is.null(   names(          qperm(app_qdr_m,      p =       0.75    ))                                            ))
stopifnot(identical(   names(          qperm(app_qdr_m,      p = c(s = 0.75)   )), "s"                                       )) # unnamed in < 1.3-0
stopifnot(  is.null(   names(          rperm(app_qdr_m,      n =       5       ))                                            ))
stopifnot(  is.null(   names(         pvalue(app_qdr_m, method = "global"      ))                                            ))
stopifnot(  is.null(   names(         pvalue(app_qdr_m, method = "single-step" ))                                            ))
stopifnot(  is.null(   names(         pvalue(app_qdr_m, method = "step-down"   ))                                            ))
stopifnot(  is.null(   names(         pvalue(app_qdr_m, method = "unadjusted"  ))                                            ))
stopifnot(  is.null(   names(      midpvalue(app_qdr_m                         ))                                            ))
stopifnot(identical(   names(pvalue_interval(app_qdr_m                         )), c("p_0", "p_1")                           ))
stopifnot(  is.null(   names(           size(app_qdr_m,  alpha = 0.05          ))                                            ))


###
### Exact using shift algorithm
###

### Scalar
shf_scl <- independence_test(y1 ~ x, distr = exact(algo = "shift"))
s <- statistic(shf_scl)

stopifnot(identical(   names(    expectation(shf_scl                         )), "1"            ))
stopifnot(identical(   names(       variance(shf_scl                         )), "1"            ))
stopifnot(identical(dimnames(     covariance(shf_scl                         )), list("1", "1") ))
stopifnot(  is.null(   names(      statistic(shf_scl,   type = "test"        ))                 )) # named in < 1.3-0
stopifnot(identical(dimnames(      statistic(shf_scl,   type = "linear"      )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(shf_scl,   type = "centered"    )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(shf_scl,   type = "standardized")), list("1", "")  ))
stopifnot(  is.null(   names(        support(shf_scl                         ))                 ))
stopifnot(  is.null(   names(          dperm(shf_scl,      x =       s       ))                 ))
stopifnot(identical(   names(          dperm(shf_scl,      x = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          pperm(shf_scl,      q =       s       ))                 ))
stopifnot(identical(   names(          pperm(shf_scl,      q = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          qperm(shf_scl,      p =       0.75    ))                 ))
stopifnot(identical(   names(          qperm(shf_scl,      p = c(s = 0.75)   )), "s"            ))
stopifnot(  is.null(   names(          rperm(shf_scl,      n =       5       ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_scl, method = "global"      ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_scl, method = "single-step" ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_scl, method = "step-down"   ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_scl, method = "unadjusted"  ))                 ))
stopifnot(  is.null(   names(      midpvalue(shf_scl                         ))                 ))
stopifnot(identical(   names(pvalue_interval(shf_scl                         )), c("p_0", "p_1")))
stopifnot(  is.null(   names(           size(shf_scl,  alpha = 0.05          ))                 ))

### Quadratic univarite
shf_qdr_u <- independence_test(y1 ~ x, distr = exact(algo = "shift"), teststat = "quadratic")
s <- statistic(shf_qdr_u)

stopifnot(identical(   names(    expectation(shf_qdr_u                         )), "1"            ))
stopifnot(identical(   names(       variance(shf_qdr_u                         )), "1"            ))
stopifnot(identical(dimnames(     covariance(shf_qdr_u                         )), list("1", "1") ))
stopifnot(  is.null(   names(      statistic(shf_qdr_u,   type = "test"        ))                 ))
stopifnot(identical(dimnames(      statistic(shf_qdr_u,   type = "linear"      )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(shf_qdr_u,   type = "centered"    )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(shf_qdr_u,   type = "standardized")), list("1", "")  ))
stopifnot(  is.null(   names(        support(shf_qdr_u                         ))                 ))
stopifnot(  is.null(   names(          dperm(shf_qdr_u,      x =       s       ))                 ))
stopifnot(identical(   names(          dperm(shf_qdr_u,      x = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          pperm(shf_qdr_u,      q =       s       ))                 ))
stopifnot(identical(   names(          pperm(shf_qdr_u,      q = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          qperm(shf_qdr_u,      p =       0.75    ))                 ))
stopifnot(identical(   names(          qperm(shf_qdr_u,      p = c(s = 0.75)   )), "s"            ))
stopifnot(  is.null(   names(          rperm(shf_qdr_u,      n =       5       ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_qdr_u, method = "global"      ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_qdr_u, method = "single-step" ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_qdr_u, method = "step-down"   ))                 ))
stopifnot(  is.null(   names(         pvalue(shf_qdr_u, method = "unadjusted"  ))                 ))
stopifnot(  is.null(   names(      midpvalue(shf_qdr_u                         ))                 ))
stopifnot(identical(   names(pvalue_interval(shf_qdr_u                         )), c("p_0", "p_1")))
stopifnot(  is.null(   names(           size(shf_qdr_u,  alpha = 0.05          ))                 ))


###
### Exact using split-up algorithm
###

### Scalar
spl_scl <- independence_test(y1 ~ x, distr = exact(algo = "split"))
s <- statistic(spl_scl)

stopifnot(identical(   names(    expectation(spl_scl                         )), "1"            ))
stopifnot(identical(   names(       variance(spl_scl                         )), "1"            ))
stopifnot(identical(dimnames(     covariance(spl_scl                         )), list("1", "1") ))
stopifnot(  is.null(   names(      statistic(spl_scl,   type = "test"        ))                 )) # named in < 1.3-0
stopifnot(identical(dimnames(      statistic(spl_scl,   type = "linear"      )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(spl_scl,   type = "centered"    )), list("1", "")  ))
stopifnot(identical(dimnames(      statistic(spl_scl,   type = "standardized")), list("1", "")  ))
## stopifnot(  is.null(   names(        support(spl_scl                         ))                 )) # NA
## stopifnot(  is.null(   names(          dperm(spl_scl,      x =       s       ))                 )) # NA
## stopifnot(identical(   names(          dperm(spl_scl,      x = c(s = s)      )), "s"            )) # NA
stopifnot(  is.null(   names(          pperm(spl_scl,      q =       s       ))                 ))
stopifnot(identical(   names(          pperm(spl_scl,      q = c(s = s)      )), "s"            ))
stopifnot(  is.null(   names(          qperm(spl_scl,      p =       0.75    ))                 ))
stopifnot(identical(   names(          qperm(spl_scl,      p = c(s = 0.75)   )), "s"            ))
stopifnot(  is.null(   names(          rperm(spl_scl,      n =       5       ))                 ))
stopifnot(  is.null(   names(         pvalue(spl_scl, method = "global"      ))                 ))
stopifnot(  is.null(   names(         pvalue(spl_scl, method = "single-step" ))                 ))
stopifnot(  is.null(   names(         pvalue(spl_scl, method = "step-down"   ))                 ))
stopifnot(  is.null(   names(         pvalue(spl_scl, method = "unadjusted"  ))                 ))
## stopifnot(  is.null(   names(      midpvalue(spl_scl                         ))                 )) # NA
## stopifnot(identical(   names(pvalue_interval(spl_scl                         )), c("p_0", "p_1"))) # NA
## stopifnot(  is.null(   names(           size(spl_scl,  alpha = 0.05          ))                 )) # NA
