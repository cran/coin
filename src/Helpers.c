/**
    Some additional functionality for package 'coin'
    *\file Helpers.c
    *\author $Author: hnilsson $
    *\date $Date: 2016-11-23 21:09:24 +0100 (Mit, 23 Nov 2016) $
*/

#include "coin_common.h"

int nrow(SEXP x) {
    SEXP a;

    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) {
        return(LENGTH(x));
    } else {
        return(INTEGER(getAttrib(x, R_DimSymbol))[0]);
    }
}

int ncol(SEXP x) {
    SEXP a;

    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) {
        return(1);
    } else {
        return(INTEGER(getAttrib(x, R_DimSymbol))[1]);
    }
}


/**
    compute a permutation of a (random subset of) 0:(m-1)
    *\param x an integer vector of length m
    *\param m integer
    *\param k integer
    *\param ans an integer vector of length k
*/

void C_SampleNoReplace(int *x, int m, int k, int *ans) {

    int i, j, n = m;

    for (i = 0; i < m; i++)
        x[i] = i;
    for (i = 0; i < k; i++) {
        j = n * unif_rand();
        ans[i] = x[j];
        x[j] = x[--n];
    }
}


SEXP R_blocksetup (SEXP block) {

    int n, nlev, nlevels, i, j, *iblock, l;
    SEXP ans, dims, indices, dummies, pindices, lindex;

    /* block _must_ be a factor! */
    n = LENGTH(block);
    iblock = INTEGER(block);
    nlevels = 1;
    for (i = 0; i < n; i++) {
        if (iblock[i] > nlevels)
            nlevels = iblock[i];
    }

    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, dims = allocVector(INTSXP, 2));
    SET_VECTOR_ELT(ans, 1, indices = allocVector(VECSXP, nlevels));
    SET_VECTOR_ELT(ans, 2, dummies = allocVector(VECSXP, nlevels));
    SET_VECTOR_ELT(ans, 3, pindices = allocVector(VECSXP, nlevels));

    INTEGER(dims)[0] = n;
    INTEGER(dims)[1] = nlevels;

    for (l = 1; l <= nlevels; l++) {

        /* number of elements in block 'l' */
        nlev = 0;
        for (i = 0; i < n; i++) {
            if (iblock[i] == l) nlev++;
        }

        /* which(block == l) and memory setup */
        SET_VECTOR_ELT(indices, l - 1, lindex = allocVector(INTSXP, nlev));
        SET_VECTOR_ELT(dummies, l - 1, allocVector(INTSXP, nlev));
        SET_VECTOR_ELT(pindices, l - 1, allocVector(INTSXP, nlev));

        j = 0;
        for (i = 0; i < n; i++) {
            if (iblock[i] == l) {
                INTEGER(lindex)[j] = i;
                j++;
            }
        }
    }

    UNPROTECT(1);
    return(ans);
}


/**
    Block permutation
    *\param blocksetup as computed by 'R_blocksetup'
    *\param ans integer vector
*/

void C_blockperm (SEXP blocksetup, int *ans) {

    int nlevels, l, nlev, j, *iindex, *ipindex;
    SEXP indices, dummies, pindices, index, dummy, pindex;

    /* n = INTEGER(VECTOR_ELT(blocksetup, 0))[0]; not used*/
    nlevels = INTEGER(VECTOR_ELT(blocksetup, 0))[1];
    indices = VECTOR_ELT(blocksetup, 1);
    dummies = VECTOR_ELT(blocksetup, 2);
    pindices = VECTOR_ELT(blocksetup, 3);

    for (l = 1; l <= nlevels; l++) {

        /* number of elements in block 'l' */
        index = VECTOR_ELT(indices, l - 1);
        dummy = VECTOR_ELT(dummies, l - 1);
        pindex = VECTOR_ELT(pindices, l - 1);
        nlev = LENGTH(index);
        iindex = INTEGER(index);
        ipindex = INTEGER(pindex);

        C_SampleNoReplace(INTEGER(dummy), nlev, nlev, ipindex);

        for (j = 0; j < nlev; j++) {
            ans[iindex[j]] = iindex[ipindex[j]];
        }
    }
}

SEXP R_blockperm (SEXP block) {

    SEXP blocksetup, ans;

    blocksetup = R_blocksetup(block);
    PROTECT(ans = allocVector(INTSXP, LENGTH(block)));
    GetRNGstate();
    C_blockperm(blocksetup, INTEGER(ans));
    PutRNGstate();
    UNPROTECT(1);
    return(ans);
}

SEXP R_MonteCarloIndependenceTest (SEXP x, SEXP y, SEXP block, SEXP B) {

    int n, p, q, pq, i, *index, *permindex, b, Bsim;
    SEXP ans, blocksetup, linstat;
    double *dans, *dlinstat, *dx, *dy, f = 0.1;

    n = nrow(x);
    p = ncol(x);
    q = ncol(y);
    pq = p*q;
    Bsim = INTEGER(B)[0];
    dx = REAL(x);
    dy = REAL(y);

    index = Calloc(n, int);
    permindex = Calloc(n, int);

    PROTECT(blocksetup = R_blocksetup(block));

    PROTECT(ans = allocMatrix(REALSXP, pq, Bsim));
    dans = REAL(ans);
    PROTECT(linstat = allocVector(REALSXP, pq));
    dlinstat = REAL(linstat);

    for (i = 0; i < n; i++)
        index[i] = i;

    GetRNGstate();

    for (b = 0; b < Bsim; b++) {

        C_blockperm(blocksetup, permindex);
        C_PermutedLinearStatistic(dx, p, dy, q, n, n, index, permindex, dlinstat);

        for (i = 0; i < pq; i++) dans[b*pq + i] = dlinstat[i];

        /* check user interrupts */
        if (b > Bsim * f) {
            R_CheckUserInterrupt();
            f += 0.1;
        }
    }

    PutRNGstate();

    Free(index); Free(permindex);

    UNPROTECT(3);
    return(ans);
}


SEXP R_maxstattrafo(SEXP x, SEXP cutpoints) {

    int i, j, n, nc, jn;
    SEXP ans;
    double *dans, *dx, *dcutpoints, cj;

    if (!isReal(x) || !isReal(cutpoints))
        error("x or cutpoints are not of type REALSXP");

    n = LENGTH(x);
    nc = LENGTH(cutpoints);
    PROTECT(ans = allocMatrix(REALSXP, n, nc));
    dans = REAL(ans);
    dx = REAL(x);
    dcutpoints = REAL(cutpoints);

    for (j = 0; j < nc; j++) {
        jn = j * n;
        cj = dcutpoints[j];
        for (i = 0; i < n; i++) {
            if (dx[i] > cj) {
                dans[jn + i] = 0.0;
            } else {
                dans[jn + i] = 1.0;
            }
        }
    }
    UNPROTECT(1);
    return(ans);
}
