/**
    Some additional functionality for package 'coin'

    *\file Helpers.c
*/

#include "coin_common.h"
#include <R_ext/Rdynload.h>
#include <libcoinAPI.h>

SEXP R_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset,
                                      SEXP block, SEXP varonly, SEXP tol) {
    return(libcoin_R_ExpectationCovarianceStatistic(x, y, weights, subset,
                                                    block, varonly, tol));
}

SEXP R_PermutedLinearStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset,
                               SEXP block, SEXP nresample) {
    return(libcoin_R_PermutedLinearStatistic(x, y, weights, subset,
                                             block, nresample));
}

SEXP R_quadform(SEXP linstat, SEXP expect, SEXP MPinv_sym) {
    return(libcoin_R_quadform(linstat, expect, MPinv_sym));
}

SEXP R_kronecker(SEXP A, SEXP B) {
    return(libcoin_R_kronecker(A, B));
}

SEXP R_MPinv_sym(SEXP x, SEXP n, SEXP tol) {
    return(libcoin_R_MPinv_sym(x, n, tol));
}

SEXP R_unpack_sym(SEXP x, SEXP names, SEXP diagonly) {
    return(libcoin_R_unpack_sym(x, names, diagonly));
}

int nrow(SEXP x) {
    SEXP a;

    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) {
        return(LENGTH(x));
    } else {
        return(INTEGER(a)[0]);
    }
}

int ncol(SEXP x) {
    SEXP a;

    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) {
        return(1);
    } else {
        return(INTEGER(a)[1]);
    }
}

SEXP R_maxstattrafo(SEXP x, SEXP cutpoints) {

    int n, nc, jn;
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

    for (int j = 0; j < nc; j++) {
        jn = j * n;
        cj = dcutpoints[j];
        for (int i = 0; i < n; i++) {
            if (ISNAN(dx[i])) {
                dans[jn + i] = dx[i];
            } else if (dx[i] > cj) {
                dans[jn + i] = 0.0;
            } else {
                dans[jn + i] = 1.0;
            }
        }
    }
    UNPROTECT(1);
    return(ans);
}

/**
    Computes the outer sum of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_outersum (const double *A, const int m, const int n,
                 const double *B, const int r, const int s,
                 double *ans) {

    int mr, js, ir;
    double y;

    mr = m * r;
    for (int i = 0; i < m; i++) {
        ir = i * r;
        for (int j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (int k = 0; k < r; k++) {
                for (int l = 0; l < s; l++) {
                    ans[(js + l) * mr + ir + k] = y + B[l * r + k];
                }
            }
        }
    }
}


/**
    R-interface to C_outersum\n
    *\param A matrix
    *\param B matrix
*/

SEXP R_outersum(SEXP A, SEXP B) {

    int m, n, r, s;
    SEXP ans;

    if (!isReal(A) || !isReal(B))
        error("R_outersum: A and / or B are not of type REALSXP");

    m = nrow(A);
    n = ncol(A);
    r = nrow(B);
    s = ncol(B);

    PROTECT(ans = allocVector(REALSXP, m * n * r * s));
    C_outersum(REAL(A), m, n, REAL(B), r, s, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
