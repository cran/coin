
/**
    Exact Distribution of Two-Sample Permutation Tests   
    Streitberg & Roehmel Algorithm

    *\file StreitbergRoehmel.c
    *\author $Author: hothorn $
    *\date $Date: 2005/07/28 15:04:29 $
*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

/*
	length(scores) <= 1.000.000 observations only.
*/

#define PERM_MAX_N 1000000


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
      (usualy score_a = rep(1, length(score_a)) and 
              score_b = Data scores, Wilcoxon, Ansari ...).
      In this case the exact conditional distribution 
      in the simple independent two-sample problem is computed.
    */ 

    int n, im_a, im_b;		/* number of observations */

    SEXP H, x;		        /* matrix of permutations and vector 
                                   of probabilities */ 
  
    int i, j, k, sum_a = 0, sum_b = 0, s_a = 0, s_b = 0, isb;
    double msum = 0.0; 	        /* little helpers */
  
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

    if (n > PERM_MAX_N)
        error("n > %d in R_cpermdistr2", PERM_MAX_N); 

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
            dx[j] = dH[isb + j + 1];
            msum += dx[j];
        }
	
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
