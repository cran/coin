/**
    S4 classes
    *\file Classes.c
    *\author $Author: hnilsson $
    *\date $Date: 2015-07-27 18:00:00 +0200 (Mon, 27 Jul 2015) $
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
