/**
    Exact Distribution of Two-Sample Permutation Tests
    Streitberg & Roehmel Algorithm

    *\file StreitbergRoehmel.c
    *\author $Author: hnilsson $
    *\date $Date: 2015-07-27 18:00:00 +0200 (Mon, 27 Jul 2015) $
*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

/**
    The density of the permutation distribution for
    the independent two sample problem.

    REFERENCES

    Bernd Streitberg & Joachim R\"ohmel (1986),
    Exact distributions for permutations and rank tests:
    An introduction to some recently published algorithms.
    Statistical Software Newsletter 12(1), 10-17.

    Bernd Streitberg & Joachim R\"ohmel (1987),
    Exakte Verteilungen f\"ur Rang- und Randomisierungstests
    im allgemeinen $c$-Stichprobenfall.
    EDV in Medizin und Biologie 18(1), 12-19 (in german).

    *\param score_a score vector (typically c(1,1,...,1))
    *\param score_b score vector (typically ranks)
    *\param m_a integer indicating the sum of m_a elements of score_a
    *\param m_b integer indicating the sum of m_b elements of score_b
    *\param retProb logical indicating whether the density (TRUE) or
            the matrix of all permutations should be returned
*/

SEXP R_cpermdist2(SEXP score_a, SEXP score_b, SEXP m_a,  SEXP m_b,
                  SEXP retProb) {
    /*
      compute the joint permutation distribution of the
      sum of the first m_a elements of score_a and score_b
      (usually score_a = rep(1, length(score_a)) and
               score_b = Data scores, Wilcoxon, Ansari ...).
      In this case the exact conditional distribution
      in the simple independent two-sample problem is computed.
    */

    int n, im_a, im_b;          /* number of observations */

    SEXP H, x;                  /* matrix of permutations and vector
                                   of probabilities */

    int i, j, k, sum_a = 0, sum_b = 0, s_a = 0, s_b = 0, isb;
    double msum = 0.0;          /* little helpers */

    int *iscore_a, *iscore_b;   /* pointers to R structures */
    double *dH, *dx;

    /* some basic checks, should be improved */

    if (!isVector(score_a))
        error("score_a is not a vector");

    n = LENGTH(score_a);

    if (!isVector(score_b))
        error("score_b is not a vector");

    if (LENGTH(score_b) != n)
        error("length of score_a and score_b differ");

    iscore_a = INTEGER(score_a);
    iscore_b = INTEGER(score_b);

    if (TYPEOF(retProb) != LGLSXP)
        error("retProb is not a logical");

    im_a = INTEGER(m_a)[0];  /* cosmetics only */
    im_b = INTEGER(m_b)[0];

    /* compute the total sum of the scores and check if they are >= 0 */

    for (i = 0; i < n; i++) {
        if (iscore_a[i] < 0)
            error("score_a for observation number %d is negative", i);
        if (iscore_b[i] < 0)
            error("score_b for observation number %d is negative", i);
        sum_a += iscore_a[i];
        sum_b += iscore_b[i];
    }

    /*
      optimization according to Streitberg & Roehmel
    */

    sum_a = imin2(sum_a, im_a);
    sum_b = imin2(sum_b, im_b);

    /*
        initialize H
    */

    PROTECT(H = allocVector(REALSXP, (sum_a + 1) * (sum_b + 1)));
    dH = REAL(H);

    for (i = 0; i <= sum_a; i++) {
        isb = i * (sum_b + 1);
        for (j = 0; j <= sum_b; j++) dH[isb + j] = 0.0;
    }

    /*
        start the Shift-Algorithm with H[0,0] = 1
    */

    dH[0] = 1.0;

    for (k = 0; k < n; k++) {
        s_a += iscore_a[k];
        s_b += iscore_b[k];

        /*
            compute H up to row im_aand column im_b
            Note:
            sum_a = min(sum_a, m)
            sum_b = min(sum_b, c)
        */

        for (i = imin2(im_a, s_a); i >= iscore_a[k]; i--) {
            isb = i * (sum_b + 1);
            for (j = imin2(im_b,s_b); j >= iscore_b[k]; j--)
                dH[isb + j] +=
                    dH[(i - iscore_a[k]) * (sum_b + 1) + (j - iscore_b[k])];
        }
    }

    /*
        return the whole matrix H
        Note: use matrix(H, nrow=m_a+1, byrow=TRUE) in R
    */

    if (!LOGICAL(retProb)[0]) {
        UNPROTECT(1);
        return(H);
    } else {
        PROTECT(x = allocVector(REALSXP, sum_b));
        dx = REAL(x);

        /*
            get the values for sample size im_a (in row m) and sum it up
        */

        isb = im_a * (sum_b + 1);
        for (j = 0; j < sum_b; j++) {
            if (!R_FINITE(dH[isb + j + 1]))
                error("overflow error; cannot compute exact distribution");
            dx[j] = dH[isb + j + 1];
            msum += dx[j];
        }
        if (!R_FINITE(msum) || msum == 0.0)
            error("overflow error; cannot compute exact distribution");

        /*
            compute probabilities and return the density x to R
            the support is min(score_b):sum(score_b)
            [dpq] stuff is done in R
        */

        for (j = 0; j < sum_b; j++)
            dx[j] = dx[j]/msum;

        UNPROTECT(2);
        return(x);
    }
}

/**
    The density of the permutation distribution for
    the one sample problem.

    REFERENCES

    Bernd Streitberg & Joachim R\"ohmel (1986),
    Exact distributions for permutations and rank tests:
    An introduction to some recently published algorithms.
    Statistical Software Newsletter 12(1), 10-17.

    Bernd Streitberg & Joachim R\"ohmel (1987),
    Exakte Verteilungen f\"ur Rang- und Randomisierungstests
    im allgemeinen $c$-Stichprobenfall.
    EDV in Medizin und Biologie 18(1), 12-19 (in german).

    *\param scores score vector (such as rank(abs(y)) for wilcoxsign_test)
*/


SEXP R_cpermdist1(SEXP scores) {

    /*
      compute the permutation distribution of the sum of the
      absolute values of the positive elements of `scores'
    */

    int n;      /* number of observations */
    SEXP H;     /* vector giving the density of statistics 0:sum(scores) */

    int i, k, sum_a = 0, s_a = 0; /* little helpers */
    int *iscores;
    double msum = 0.0;
    double *dH;

    n = LENGTH(scores);
    iscores = INTEGER(scores);

    for (i = 0; i < n; i++) sum_a += iscores[i];

    /*
      Initialize H
    */

    PROTECT(H = allocVector(REALSXP, sum_a + 1));
    dH = REAL(H);
    for (i = 0; i <= sum_a; i++) dH[i] = 0.0;

    /*
      start the shift-algorithm with H[0] = 1.0
    */

    dH[0] = 1.0;

    for (k = 0; k < n; k++) {
        s_a = s_a + iscores[k];
            for (i = s_a; i >= iscores[k]; i--)
                dH[i] = dH[i] + dH[i - iscores[k]];
    }


    /*
        get the number of permutations
    */

    for (i = 0; i <= sum_a; i++) {
        if (!R_FINITE(dH[i]))
            error("overflow error: cannot compute exact distribution");
        msum += dH[i];
    }
    if (!R_FINITE(msum) || msum == 0.0)
        error("overflow error: cannot compute exact distribution");

    /*
        compute probabilities and return the density H to R
        [dpq] stuff is done in R
    */

    for (i = 0; i <= sum_a; i++)
        dH[i] = dH[i]/msum;     /* 0 is a possible realization */

    UNPROTECT(1);
    return(H);
}
