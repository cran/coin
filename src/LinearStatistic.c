
/**
    Linear statistics for conditional inference
    *\file $RCSfile: LinearStatistic.c,v $
    *\author $Author: hothorn $
    *\date $Date: 2005/02/10 14:54:09 $
*/
    
#include "CI_common.h"

/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans) {

    int i, j, k, l, mr, js, ir;
    double y;

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++) {
                    ans[(js + l) * mr + ir + k] = y * B[l * r + k];
                }
            }
        }
    }
}  


/**
    Conditional expectation and covariance of the influence function\n
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; an object of class `ExpectCovarInfluence'
*/

void C_ExpectCovarInfluence(const double* y, const int q,
                            const double* weights, const int n, 
                            SEXP ans) {

    int i, j, k, jq;
    
    /* pointers to the slots of object ans */
    double *dExp_y, *dCov_y, *dsweights, tmp;
    
    /*  return values: set to zero initially */
    dExp_y = REAL(GET_SLOT(ans, CI_expectationSym));
    for (j = 0; j < q; j++) dExp_y[j] = 0.0;
    
    dCov_y = REAL(GET_SLOT(ans, CI_covarianceSym));
    for (j = 0; j < q*q; j++) dCov_y[j] = 0.0;
    
    dsweights = REAL(GET_SLOT(ans, CI_sumweightsSym));

    /*  compute the sum of the weights */
        
    dsweights[0] = 0;
    for (i = 0; i < n; i++) dsweights[0] += weights[i];
    if (dsweights[0] <= 1) 
        error("C_ExpectCovarInfluence: sum of weights is less than one");

    /*
     *    Expectation of the influence function
     */

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute */
    
        if (weights[i] == 0.0) continue;
    
        for (j = 0; j < q; j++)
            dExp_y[j] += weights[i] * y[j * n + i];
    }

    for (j = 0; j < q; j++)
        dExp_y[j] = dExp_y[j] / dsweights[0];


    /*
     *    Covariance of the influence function
     */

    for (i = 0; i < n; i++) {

        if (weights[i] == 0.0) continue;
     
        for (j = 0; j < q; j++) {
            tmp = weights[i] * (y[j * n + i] - dExp_y[j]);
            jq = j * q;
            for (k = 0; k < q; k++)
                dCov_y[jq + k] += tmp * (y[k * n + i] - dExp_y[k]);
        }
    }

    for (j = 0; j < q*q; j++)
        dCov_y[j] = dCov_y[j] / dsweights[0];
}


/**
    R-interface to C_ExpectCovarInfluence\n
    *\param y values of the influence function
    *\param weights case weights
*/

SEXP R_ExpectCovarInfluence(SEXP y, SEXP weights) {

    SEXP ans;
    int q, n;
    
    if (!isReal(y) || !isReal(weights))
        error("R_ExpectCovarInfluence: arguments are not of type REALSXP");
    
    n = nrow(y);
    q = ncol(y);
    
    if (LENGTH(weights) != n) 
        error("R_ExpectCovarInfluence: vector of case weights does not have %d elements", n);

    /*  allocate storage for return values */
    PROTECT(ans = NEW_OBJECT(MAKE_CLASS("ExpectCovarInfluence")));
    SET_SLOT(ans, CI_expectationSym, 
             PROTECT(allocVector(REALSXP, q)));
    SET_SLOT(ans, CI_covarianceSym, 
             PROTECT(allocMatrix(REALSXP, q, q)));
    SET_SLOT(ans, CI_sumweightsSym, 
             PROTECT(allocVector(REALSXP, 1)));

    C_ExpectCovarInfluence(REAL(y), q, REAL(weights), n, ans);
    
    UNPROTECT(4);
    return(ans);
}


/**
    Conditional expectation and covariance of the a linear statistic\n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param expcovinf an object of class `ExpectCovarInfluence'
    *\param ans return value; an object of class `ExpectCovar'
*/

void C_ExpectCovarLinearStatistic(const double* x, const int p, 
                                  const double* y, const int q,
                                  const double* weights, const int n,
                                  const SEXP expcovinf, SEXP ans) {

    int i, j, k, pq, ip;
    double sweights = 0.0, f1, f2, tmp;
    double *swx, *CT1, *CT2, *Covy_x_swx, 
           *dExp_y, *dCov_y, *dExp_T, *dCov_T;
    
    pq   = p * q;
    
    /* the expectation and covariance of the influence function */
    dExp_y = REAL(GET_SLOT(expcovinf, CI_expectationSym));
    dCov_y = REAL(GET_SLOT(expcovinf, CI_covarianceSym));
    sweights = REAL(GET_SLOT(expcovinf, CI_sumweightsSym))[0];

    if (sweights <= 1.0) 
        error("C_ExpectCovarLinearStatistic: sum of weights is less than one");

    /* prepare for storing the results */
    dExp_T = REAL(GET_SLOT(ans, CI_expectationSym));
    dCov_T = REAL(GET_SLOT(ans, CI_covarianceSym));

    /* allocate storage: all helpers, initially zero */
    swx = Calloc(p, double);               /* p x 1  */
    CT1 = Calloc(p * p, double);           /* p x p  */

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute */
        if (weights[i] == 0.0) continue;
    
        ip = i*p;
        for (k = 0; k < p; k++) {
            tmp = weights[i] * x[k * n + i];
            swx[k] += tmp;

            /* covariance part */
            for (j = 0; j < p; j++) {
                CT1[j * p + k] += tmp * x[j * n + i];
            }
        }
    }

    /*
    *   dExp_T: expectation of the linear statistic T
    */

    for (k = 0; k < p; k++) {
        for (j = 0; j < q; j++)
            dExp_T[j * p + k] = swx[k] * dExp_y[j];
    }

    /* 
    *   dCov_T:  covariance of the linear statistic T
    */

    f1 = sweights/(sweights - 1);
    f2 = (1/(sweights - 1));

    if (pq == 1) {
        dCov_T[0] = f1 * dCov_y[0] * CT1[0];
        dCov_T[0] -= f2 * dCov_y[0] * swx[0] * swx[0];
    } else {
        /* two more helpers needed */
        CT2 = Calloc(pq * pq, double);            /* pq x pq */
        Covy_x_swx = Calloc(pq * q, double);      /* pq x q  */
        
        C_kronecker(dCov_y, q, q, CT1, p, p, dCov_T);
        C_kronecker(dCov_y, q, q, swx, p, 1, Covy_x_swx);
        C_kronecker(Covy_x_swx, pq, q, swx, 1, p, CT2);

        for (k = 0; k < (pq * pq); k++)
            dCov_T[k] = f1 * dCov_T[k] - f2 * CT2[k];

        /* clean up */
        Free(CT2); Free(Covy_x_swx);
    }

    /* clean up */
    Free(swx); Free(CT1); 
}


/**
    R-interface to C_ExpectCovarLinearStatistic\n
    *\param x values of the transformation
    *\param y values of the influence function
    *\param weights case weights
    *\param expcovinf an object of class `ExpectCovarInfluence'
*/

SEXP R_ExpectCovarLinearStatistic(SEXP x, SEXP y, SEXP weights, 
                                  SEXP expcovinf) {
    
    SEXP ans;
    int n, p, q, pq;

    /* determine the dimensions and some checks */

    n  = nrow(x);
    p  = ncol(x);
    q  = ncol(y);
    pq = p * q;
    
    if (nrow(y) != n)
        error("y does not have %d rows", n);
    if (LENGTH(weights) != n) 
        error("vector of case weights does not have %d elements", n);

    PROTECT(ans = NEW_OBJECT(MAKE_CLASS("ExpectCovar")));
    SET_SLOT(ans, CI_expectationSym, 
             PROTECT(allocVector(REALSXP, pq)));
    SET_SLOT(ans, CI_covarianceSym, 
             PROTECT(allocMatrix(REALSXP, pq, pq)));

    C_ExpectCovarLinearStatistic(REAL(x), p, REAL(y), q, 
        REAL(weights), n, expcovinf, ans);
    
    UNPROTECT(3);
    return(ans);
}

/**
    Computes the linear statistic, formula (1) in the paper\n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/
  
void C_LinearStatistic (const double *x, const int p,
                        const double *y, const int q,
                        const double *weights, const int n,
                        double *ans) {
              
    int i, j, k, kp, kn, ip;
    double tmp;

    for (k = 0; k < q; k++) {

        kn = k * n;
        kp = k * p;
        for (j = 0; j < p; j++) ans[kp + j] = 0.0;
            
        for (i = 0; i < n; i++) {
                
            /* optimization: weights are often zero */
            if (weights[i] == 0.0) continue;
                
            tmp = y[kn + i] * weights[i];
                
            ip = i * p;
            for (j = 0; j < p; j++)
                 ans[kp + j] += x[j*n + i] * tmp;
        }
    }
}


/**
    R-interface to C_LinearStatistic \n
    *\param x values of the transformation
    *\param y values of the influence function
    *\param weights case weights
*/

SEXP R_LinearStatistic(SEXP x, SEXP y, SEXP weights) {

    /* the return value; a vector of type REALSXP */
    SEXP ans;

    /* dimensions */
    int n, p, q;

    /* 
     *    only a basic check: we do not coerce objects since this
     *    function is for internal use only
     */
    
    if (!isReal(x) || !isReal(y) || !isReal(weights))
        error("LinStat: arguments are not of type REALSXP");
    
    n = nrow(y);
    if (nrow(x) != n || LENGTH(weights) != n)
        error("LinStat: dimensions don't match");

    q    = ncol(y);
    p    = ncol(x);
           
    PROTECT(ans = allocVector(REALSXP, p*q));
 
    C_LinearStatistic(REAL(x), p, REAL(y), q, REAL(weights), n, 
                      REAL(ans));

    UNPROTECT(1);
    return(ans);
}


/**
    Linear Statistic with permuted indices\n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param n number of observations
    *\param nperm number of permutations
    *\param indx indices for the x-part
    *\param perm (permuted) indices for the y-part
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/

void C_PermutedLinearStatistic(const double *x, const int p,
                               const double *y, const int q,
                               const int n, const int nperm,
                               const int *indx, const int *perm, 
                               double *ans) {

    int i, j, k, kp, kn, knpi;

    for (k = 0; k < q; k++) {

        kn = k * n;
        kp = k * p;
        for (j = 0; j < p; j++) ans[kp + j] = 0.0;
            
        for (i = 0; i < nperm; i++) {
                
            knpi = kn + perm[i];

            for (j = 0; j < p; j++)
                ans[kp + j] += x[j*n + indx[i]] * y[knpi];
        }
    }
}


/**
    Linear Statistic with permuted indices\n
    *\param x values of the transformation
    *\param y values of the influence function
    *\param indx indices for the x-part
    *\param perm (permuted) indices for the y-part
*/

SEXP R_PermutedLinearStatistic(SEXP x, SEXP y, SEXP indx, SEXP perm) {

    SEXP ans;
    int n, nperm, p, q, i, *iperm, *iindx;

    /* 
       only a basic check
    */

    if (!isReal(x) || !isReal(y))
        error("R_PermutedLinearStatistic: arguments are not of type REALSXP");
    
    if (!isInteger(perm))
        error("R_PermutedLinearStatistic: perm is not of type INTSXP");
    if (!isInteger(indx))
        error("R_PermutedLinearStatistic: indx is not of type INTSXP");
    
    n = nrow(y);
    nperm = LENGTH(perm);
    iperm = INTEGER(perm);
    if (LENGTH(indx)  != nperm)
        error("R_PermutedLinearStatistic: dimensions don't match");
    iindx = INTEGER(indx);

    if (nrow(x) != n)
        error("R_PermutedLinearStatistic: dimensions don't match");

    for (i = 0; i < nperm; i++) {
        if (iperm[i] < 0 || iperm[i] > (n - 1) )
            error("R_PermutedLinearStatistic: perm is not between 1 and nobs");
        if (iindx[i] < 0 || iindx[i] > (n - 1) )
            error("R_PermutedLinearStatistic: indx is not between 1 and nobs");
    }

    q    = ncol(y);
    p    = ncol(x);
           
    PROTECT(ans = allocVector(REALSXP, p*q));
    
    C_PermutedLinearStatistic(REAL(x), p, REAL(y), q, n, nperm,
                 iindx, iperm, REAL(ans));
    
    UNPROTECT(1);
    return(ans);
}


/**
    Score matrix for a linear combination of a linear statistic 
    from the left\n
    *\param x score vector of length p
    *\param p length of x
    *\param q dimension of the influence function
    *\param ans return value; a pointer to a q x (pq) matrix 
*/ 

void C_scmatleft(const double *x, const int p, 
                 const int q, double *ans) {

    /*
     *
     *    The basic difficulty is that the statistic itself is
     *    L ~ (p, q) but we only look at vec(L)
     *    
     *    The function creates a score matrix A in order to
     *    compute A %*% vec(L) == x %*% L
     *
     *    pq is the length of vec(L) (and p is determined from 
     *    the length of x)
     */
    
    int k, j, pq;
    
    pq = p * q;
    for (j = 0; j < q; j++) {
            for (k = 0; k < p; k++) {
                ans[pq * j + q*k +  j] = x[k];
            }
    }
}


/**
    R-interface to C_scmatleft
    *\param x score vector of length p
    *\param pq dimension of the linear statistic
*/ 

SEXP R_scmatleft(SEXP x, SEXP pq) {

    SEXP ans;
    double *dans, *dx;
    int p, q, i;
    
    if (!isReal(x)) error("R_scmatleft: x not of type REALSXP");
    if (!isInteger(pq)) error("R_scmatleft: pq not of type INTSXP");
    
    dx = REAL(x);
    p = LENGTH(x);
    q = INTEGER(pq)[0] / p;
    
    PROTECT(ans = allocMatrix(REALSXP, q, p*q));
    dans = REAL(ans);
    for (i = 0; i < q*p*q; i++) dans[i] = 0.0;
    
    C_scmatleft(dx, p, q, dans);
    
    UNPROTECT(1);
    return(ans);
}


/**
    Score matrix for a linear combination of a linear statistic 
    from the right\n
    *\param x score vector of length q
    *\param p dimension of the transformation
    *\param q length of x 
    *\param ans return value; a pointer to a p x (pq) matrix 
*/ 

void C_scmatright(const double *x, const int p, 
                 const int q, double *ans) { 

    /* 
     *
     *    The basic difficulty is that the statistic itself is
     *    L ~ (p, q) but we only look at vec(L)
     *    
     *    The function creates a score matrix A in order to
     *    compute A %*% vec(L) == L %*% x
     *
     *    pq is the length of vec(L) (and q is determined from 
     *    the length of x)
     *
     */ 

    int i, k, pp;
    
    pp = p * p;
    for (k = 0; k < q; k++) {
        for (i = 0; i < p; i++) {
            ans[pp * k + i * p   + i] = x[k];
        }
    }
}

/**
    R-interface to C_scmatright
    *\param x score vector of length q
    *\param pq dimension of the linear statistic
*/ 

SEXP R_scmatright(SEXP x, SEXP pq) {

    SEXP ans;
    double *dans, *dx;
    int p, q, i;
    
    if (!isReal(x)) error("R_scmatright: x not of type REALSXP");
    if (!isInteger(pq)) error("R_scmatright: pq not of type INTSXP");
    
    dx = REAL(x);
    q = LENGTH(x);
    p = INTEGER(pq)[0] / q;
    
    PROTECT(ans = allocMatrix(REALSXP, p, p*q));
    dans = REAL(ans);
    
    for (i = 0; i < p*p*q; i++) dans[i] = 0.0;
    
    C_scmatright(dx, p, q, dans);

    UNPROTECT(1);
    return(ans);
}
