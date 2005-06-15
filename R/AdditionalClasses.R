
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

