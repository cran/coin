
/**
    S4 classes 
    *\file $RCSfile: Classes.c,v $
    *\author $Author: hothorn $
    *\date $Date: 2005/02/10 14:54:09 $
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
