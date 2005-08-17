
/**
    S4 classes 
    *\file Classes.c
    *\author $Author: hothorn $
    *\date $Date: 2005/07/28 15:04:29 $
*/

#include "CI_common.h"

SEXP 
    CI_expectationSym,
    CI_covarianceSym,
    CI_sumweightsSym;

SEXP coin_init(void) {
    CI_expectationSym = install("expectation");
    CI_covarianceSym = install("covariance");
    CI_sumweightsSym = install("sumweights");
    return(R_NilValue);
}
