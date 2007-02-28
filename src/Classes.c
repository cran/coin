
/**
    S4 classes 
    *\file Classes.c
    *\author $Author: hothorn $
    *\date $Date: 2007-02-15 09:25:46 +0100 (Thu, 15 Feb 2007) $
*/

#include "coin_common.h"

SEXP 
    coin_expectationSym,
    coin_covarianceSym,
    coin_sumweightsSym;

SEXP coin_init(void) {
    coin_expectationSym = install("expectation");
    coin_covarianceSym = install("covariance");
    coin_sumweightsSym = install("sumweights");
    return(R_NilValue);
}
